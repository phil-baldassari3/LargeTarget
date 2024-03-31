"""
fitness4gene_size version 1.6

Author: Phil Baldassari

Description: Wright-Fisher process on two differently sized genes. One gene has a longer CDS and more TFBS than the other.
In this simplistic model, a gene is expressed when at least one TFBS in the CRM is bound by a TF.
A longer CDS has more chance for Non-synonymous mutation.

Model:
Each gene has a certain number of TFBS each with a probability score of it being bound by a TF.
The Cis-regulatory model score cooresponds to the probability that at least one TFBS was bound at a given moment.
The CDS recieves an initial score of 1 meaning a normal gene product.
A score of 1 indicates a normal gene product with lower scores being less functional and higher scores being more functional.
Each generation mutations are added based on a per site mutation rate. The number of mutations is sampled from a binomial distribution.
Mutations in TFBS affect the TFBS score propotianally. Multiplpiers are sampled from a gamma distribution such that a mutation increases the 
binding probability 10% of the time and on average the binding affinity decreases to ~60% of the original affinity.
Mutations on CDS only affect function if they are nonsynonymous which occurs with around 73% of all mutations in CDS. These muations also affect function propotianlly.
Mutational effect proportions are sampled from a gamma distribution such that the expected effect of a NS mutation decreases the gene function
by ~65% and mutations will be at least neutral ~1% of the time.
The Cis-regulatory model scores are analogous to a relative expression level. These scores are input into one of 3 fitness functions for selection at each genreation: i) linear, ii) parabolic, iii) sigmoidal
to compute the CRM's fitness contribution score. The gene score is computed by multiplying the CRM fitness score and the CDS score.
Since species is diploid, the scores of both homologs are averaged for an individual score.
The population of individuals are diploid and obligate outbreeding hermaphrodites. The genreation size remains constant.

usage: Use through the command line. See README.md for details.
"""

#importing modules
import copy
import math
import random
import numpy as np
import pandas as pd
from statistics import mean
from statistics import variance
import matplotlib.pyplot as plt
import argparse



#parsing command line arguments
parser = argparse.ArgumentParser()

#arguments
parser.add_argument('--mu', type=float, default=0.001, help='per site mutation rate, default is 0.001')
parser.add_argument('--pop', type=int, default=100, help='constant population size, default is 100')
parser.add_argument('--gen', type=int, default=100, help='number of genreations for simulation, default is 100')
parser.add_argument('--express2fit', type=str, required=True, help='function repressenting relationship b/n expression and fitness, options are: "linear", "parabolic", or "sigmoidal"')
parser.add_argument('--burnin', type=int, default=50, help='number of burnin generations, default is 100')
parser.add_argument('--seed', type=int, required=True, help='seed value for randomization modules')
parser.add_argument('--tfbs1', type=int, required=True, help='number of TFBSs for gene 1')
parser.add_argument('--tfbslen1', type=int, required=True, help='average length of TFBSs for gene 1')
parser.add_argument('--cds1', type=int, required=True, help='length of CDS for gene 1 note: this number does NOT have to be divisible by 3')
parser.add_argument('--tfbs2', type=int, required=True, help='number of TFBSs for gene 2')
parser.add_argument('--tfbslen2', type=int, required=True, help='average length of TFBSs for gene 2')
parser.add_argument('--cds2', type=int, required=True, help='length of CDS for gene 2 note: this number does NOT have to be divisible by 3')
parser.add_argument('--mueffects', type=str, default="all", help='control parameter to control the effect of mutations, default is "all", change to "deleterious" for all mutations to be deleterious')
parser.add_argument('--mutypes', type=str, default="all", help='control parameter to control if mutations land in TFBSs, CDS, or both, default is "all", change to "CDS" for all mutations to only be in CDS or "CRM" for mutations to only be in TFBSs')
parser.add_argument('--outputdir', type=str, default="", help='path to the directory in which to save the output files, note: the directory sting must end in "/" and it can be a relative path if the output directory is a subdirectory to the one in which this script is saved.')

#parse args
args = parser.parse_args()

#assign variables
mu = args.mu
n = args.pop
g = args.gen
fitness_function = args.express2fit
burnin_g = args.burnin
seed = args.seed
tfbs1 = args.tfbs1
tfbs_len1 = args.tfbslen1
CDS_len1 = args.cds1
tfbs2 = args.tfbs2
tfbs_len2 = args.tfbslen2
CDS_len2 = args.cds2
mutation_effects = args.mueffects
mutations_types = args.mutypes
out_dir = args.outputdir



#setting mutational effects sampling interval
if mutation_effects == "all":
    mueff_interval = (float('-inf'), float('inf'))
elif mutation_effects == "deleterious":
    mueff_interval = (0, 0.9)
else:
    print("mueffects parameter incorrectly set. using all mutational effects.")
    mueff_interval = (float('-inf'), float('inf'))



#fitness function parameter, starts at -1 and will be calibrated later
a = -1

#fitness functions
def linear(x):
    """f(x) = x - a"""
    y = x - a

    return y


def parabolic(x):
    """f(x) = -(2x - a)^2 + 1"""

    y = -1 * ((2*x)-a)**2 + 1

    return y


def sigmoidal(x):
    """f(x) = 1/(1 + e^-8(2x - a))"""

    y = 1 / (1 +(math.e ** (-8*(2*x - a))))

    return y

#fitness function dictionary; to be used for selecting the fitness function to use i.e. fitness_func["function"]
fitness_func = {"linear":linear, "parabolic":parabolic, "sigmoidal":sigmoidal}



#Computational functions
def bound_at_zero(ls):
    """
    Function takes in a list of scores or fitness values and bounds the list at zero. Any negative values are turned to zero.
    Returns new list
    """

    new_ls = [0 if x < 0 else x for x in ls]

    return new_ls


def normalize(ls):
    """Normalizes a set of values in a list to scale between 0 and 1. Returns a new list"""

    sumoflist = sum(ls)

    norm = [i/sumoflist for i in ls]

    return norm


def compute_CRM_scores(ls):
    """
    Function takes the list of TFBS scores and uses them to compute the probability of at least one TFBS in the CRM being bound at least once at a given moment.
    Pr(e) = 1 - PI(1-b)
    Returns the numerical value
    """
    score = 1

    for i in ls:
        score *= (1 - i)

    final_score = (1 - score)

    return final_score


def compute_indv_fitness(pop_ls):
    """
    Function takes in a population list and computes a relative fitness score for each individual by coputing a score for each gene and averaging the score for each gene pair (diploid individual)
    Scores are coputed by the product of TFBS fitness contribution and CDS fitness contribution. The TFBS fitness contribution is computed by calulating a raw score based on binding probabilities
    using the compute_CRM_scores function. This raw score is used in the chosen expression-to-fitness function to output the TFBS fitness contribution.
    The CDS fitness contribution begins at 1 in the start of the simulation and changes propotionally by being mutlplied by mutational effects sampled from a gamma distirbution.
    The fitness scores for each gene in the homologous pair are averaged and the average is appended to the indv_w list.
    Returns a list of expression scores and a list of relative fitnesses for every indidivual in the population.
    """

    indv_expression_score = []
    indv_w = []
    
    #looping through individuals
    for indv_idx in range(len(pop_ls)):
        gene_scores = []
        gene_ws = []

        #loping through each gene per individual
        for gene_idx in range(2):

            #computing score for each gene
            TFBS_rawscore = compute_CRM_scores(pop_ls[indv_idx][gene_idx][0][0])
            TFBS_w_contrib = fitness_func[fitness_function](TFBS_rawscore)
            CDS_w_contrib = pop_ls[indv_idx][gene_idx][1][0]
            
            windv = TFBS_w_contrib * CDS_w_contrib
            
            gene_scores.append(TFBS_rawscore)
            gene_ws.append(windv)

        #averaging score for each diploid pair
        indv_expression_score.append(mean(gene_scores))
        indv_w.append(mean(gene_ws))

    return indv_w, indv_expression_score



#simulation functions
def gamma_mutational_effects(shape, scale, interval=(float('-inf'), float('inf'))):
    """
    Function for sampling from a gamma distribution to add mutational effects.
    Optional interval parameter can be used to only select deleterious mutations (or only beneficial mutations if you want).
    Returns a random sample from the gamma distriobution.
    """

    sample = -1
    while True:
        s = np.random.gamma(shape, scale)
        if interval[0] <= s <= interval[1]:
            sample = s
            break
        else:
            continue

    return sample


def burnin(num_of_tfbs, length_of_tfbs, length_of_cds):
    """
    Function to run a burnin of a set amount of generations (does not need to exceed 200).
    The average of the expression score of the final burnin population will be used to calibrate the fitness functions to set the average to a fitness of 0.5 by adjusting the a parameter.
    The final population is then returned to be used in the simulation used as the starting population.
    The burnin runs under a neutral W-F process.
    Returns starting population and sets the fitness function parameter, a.
    """

    #setting parameters
    tfbs = num_of_tfbs
    tfbs_len = length_of_tfbs
    CDS_len = length_of_cds

    #running burnin
    populationBI0, scBI0, wBI0 = starting_pop(tfbs, tfbs_len, CDS_len)
    pplnBI = populationBI0
    wBI = wBI0

    for gen in range(burnin_g):
        pplnBI, scBI, wBI = next_gen(pplnBI, wBI, selection="neutral")

    #parameterizing fitness functions
    s_bar = mean(scBI)

    if fitness_function == "linear":
        param = s_bar - 0.5

    elif fitness_function == "parabolic":
        paramls = []
        param1 = (2*s_bar) + (math.sqrt(0.5))
        param2 = (2*s_bar) - (math.sqrt(0.5))
        paramls.append(param1)
        paramls.append(param2)
        param = max(paramls)

    elif fitness_function == "sigmoidal":
        param = 2 * s_bar

    else:
        print("Fitness function set incorrectly.\n")

    #setting fitness function parameter
    global a
    a = param

    return pplnBI


def starting_pop(tfbs, tfbs_len, CDS_len):
    """
    Function randomly generates the starting population's genes based on tfbs, tfbs_len, CDS_len, and n
    Each individual is represented by a list with two lists representing their two homologous genes.
    Gene: [[[TFBS scores], [TFBS lengths]], [CDS score, CDS length]]
    Individual: [gene, gene]
    Pop: [individual, individual,...]
    Returns: population list, list of expression scores, and a list of relative fitnesses
    """

    popls = []

    #tfbs lengths
    tfbs_lens = list(np.random.gamma((tfbs_len/2), scale=2, size=tfbs))
    tfbs_lens = [int(x) for x in tfbs_lens]

    #scores for TFBS binding probabilites
    tfbs_probs = list(np.random.exponential(0.1, size=tfbs))

    #getting rid of initial probabilities at zero
    tfbs_probs = [p if p > 0 else 0.0001 for p in tfbs_probs]
    #getting rid of initial probabilities above 1
    tfbs_probs = [p if p < 1 else 0.99 for p in tfbs_probs]

    #genreating genes for each individual
    for i in range(n):
        indv = []
        #each individual is diploid
        for j in range(2):
            gene = []

            #CDS score and length
            lenofcds = CDS_len
            while True:
                if (lenofcds % 3) != 0:
                    lenofcds -= 1
                else:
                    break

            cds = [1, lenofcds]

            #append to gene
            gene.append([tfbs_probs, tfbs_lens])
            gene.append(cds)

            #appending 2 homologs to individual
            indv.append(gene)

        #place individual into population
        popls.append(indv)

    #computing raw scores for each individual
    fitnesses, expression_scores = compute_indv_fitness(popls)

    #bounding at zero
    fitnesses = bound_at_zero(fitnesses)

    return popls, expression_scores, fitnesses


def next_gen(pop_ls, ws, selection="selection"):
    """
    Fuction outputs the next generation from a previous genreation using the W-F model.
    Inputs are the population list and the list of fitness scores.
    The fitness will be used to simulate mating with a skew toward more mating instances between individuals of higher fitness (selection).
    Returns a new population list and expression scores list as well as a new fitness score list.
    """

    new_pop = []

    #selecting mating pair based on relative fitness
    for i in range(n):

        if selection == "neutral":
            pair_idx = list(np.random.choice([x for x in range(len(pop_ls))], size=2, replace=False))
        else:
            pair_idx = list(np.random.choice([x for x in range(len(pop_ls))], size=2, replace=False, p=normalize(ws)))
            
        indv1 = pop_ls[pair_idx[0]]
        indv2 = pop_ls[pair_idx[1]]

        #simulating meiosis
        offspring_indv = []
        genefrom1 = copy.deepcopy(indv1[random.randint(0, 1)])
        genefrom2 = copy.deepcopy(indv2[random.randint(0, 1)])
        offspring_indv.append(genefrom1)
        offspring_indv.append(genefrom2)

        #adding mutations
        for gene in offspring_indv:

            if mutations_types == "all" or mutations_types == "CRM":
                #looping through TFBS
                for idx in range(len(gene[0][0])):
                    #how many mutations for each TFBS
                    howmanyTFBS = np.random.binomial(gene[0][1][idx], mu)

                    #changing scores proportinally: s' = s(1+x)
                    pc_prop = 1
                    for j in range(howmanyTFBS):
                        pc_prop *= gamma_mutational_effects(3, 0.1879, interval=mueff_interval)
                    
                    if howmanyTFBS == 0:
                        continue
                    else:
                        #new score
                        newscore = gene[0][0][idx] * pc_prop
                        #setting new score
                        gene[0][0][idx] = newscore

            if mutations_types == "all" or mutations_types == "CDS":
                #CDS mutations
                #how many CDS mutations
                howmanyCDS = np.random.binomial(gene[1][1], mu)

                #are they SS or NS
                SNP_eff = []
                for k in range(howmanyCDS):
                    #if you want I can explain this, but out of the 192 possible sites in all 64 codons, 139.75 of those sites would result in a AA change if subsitiuted
                    if random.random() < 0.728:
                        SNP_eff.append("NS")
                    else:
                        SNP_eff.append("SS")

                #changing scores proportinally
                cds_pc_prop = 1
                for j in range(SNP_eff.count("NS")):
                    cds_pc_prop *= gamma_mutational_effects(3, 0.119, interval=mueff_interval)
                
                if howmanyCDS == 0:
                    continue
                else:
                    #new score
                    cds_score = gene[1][0] * cds_pc_prop

                #setting cds score
                gene[1][0] = cds_score

        #adding offspring to population
        new_pop.append(offspring_indv)

    #computing raw scores for each individual
    fitnesses, expression_scores = compute_indv_fitness(new_pop)

    #bounding at zero
    fitnesses = bound_at_zero(fitnesses)

    return new_pop, expression_scores, fitnesses


def sim_generations(population0, scores0, fitnesses0):
    """
    Function run the W-F simulation. Generates plots of average fitness and max fitness per generation.
    """

    #starting population
    ppln = population0
    w = fitnesses0

    #adding first generation to plot
    generation = [0]
    max_s = []
    avg_s = []
    var_s = []
    avg_w = []
    max_w = []
    var_w = []
    avg_Dw = [0]
    max_s.append(max(scores0))
    avg_s.append(mean(scores0))
    var_s.append(variance(scores0))
    avg_w.append(mean(fitnesses0))
    max_w.append(max(fitnesses0))
    var_w.append(variance(fitnesses0))

    #setting generation breaks for progress printing purposes
    breaks = g // 10

    #running for generations 1-g
    for gen in range(g):
        ppln, sc, w = next_gen(ppln, w)

        max_s.append(max(sc))
        avg_s.append(mean(sc))
        var_s.append(variance(sc))
        avg_w.append(mean(w))
        max_w.append(max(w))
        var_w.append(variance(w))
        generation.append(gen+1)

        #printing progress
        if gen % breaks == 0:
            print("Generation {}".format(str(gen)))

    #computing average delta fitness
    avg_Dw += list(np.diff(np.array(avg_w)))

    #computing standard deviations from variances for plotting
    sd_s = np.sqrt(var_s)
    sd_w = np.sqrt(var_w)

    #plotting
    fig, axs = plt.subplots(5, 1, figsize=(10, 15))

    axs[0].plot(generation, max_s)
    axs[0].set_title("Maximum Expression Score per Generation")
    axs[0].set_xlabel("Generation")
    axs[0].set_ylabel("Max Expression Score")

    axs[1].plot(generation, avg_s)
    axs[1].fill_between(generation, (avg_s - sd_s), (avg_s + sd_s), color='gray', alpha=0.3)
    axs[1].set_title("Average Expression Score per Generation")
    axs[1].set_xlabel("Generation")
    axs[1].set_ylabel("Avg Expression Score")


    axs[2].plot(generation, max_w)
    axs[2].set_title("Maximum Fitness per Generation")
    axs[2].set_xlabel("Generation")
    axs[2].set_ylabel("Max Relative Fitness")

    axs[3].plot(generation, avg_w)
    axs[3].fill_between(generation, (avg_w - sd_w), (avg_w + sd_w), color='gray', alpha=0.3)
    axs[3].set_title("Average Fitness per Generation")
    axs[3].set_xlabel("Generation")
    axs[3].set_ylabel("Avg Relative Fitness")

    axs[4].plot(generation, avg_Dw)
    axs[4].set_title("Average Delta Fitness per Generation")
    axs[4].set_xlabel("Generation")
    axs[4].set_ylabel("Avg Delta Fitness")

    fig.subplots_adjust(hspace=0.5)

    plt.savefig('{out}WF_plot_gene{gno}.png'.format(out=out_dir, gno=genenumber))

    #saving csv
    dictionary = {
        "generation":generation, 
        "max_expression_score": max_s, "avg_expression_score": avg_s, "variance_expression_score":var_s, 
        "max_fitness": max_w, "avg_fitness": avg_w, "variance_fitness":var_w, "avg_DELTAfitness": avg_Dw
        }
    df = pd.DataFrame(dictionary)
    df.to_csv("{out}WF_data_gene{gno}.csv".format(out=out_dir, gno=genenumber), index=False)


def run_simulator(num_of_tfbs, length_of_tfbs, length_of_cds):
    """Fucntion runs the simulator for a set of paramters"""

    print("Starting WF process for gene{}".format(genenumber))
    #seed value
    random.seed(seed)
    np.random.seed(seed)

    #setting parameters
    tfbs = num_of_tfbs
    tfbs_len = length_of_tfbs
    CDS_len = length_of_cds

    print("Running Burnin...")
    #running burnin
    population0 = burnin(tfbs, tfbs_len, CDS_len)

    #computing scores and fitnesses for each individual
    fitnesses0, expression_scores0 = compute_indv_fitness(population0)

    #bounding at zero
    w0 = bound_at_zero(fitnesses0)

    print("Simulating generations...")
    #simulating generations
    sim_generations(population0, expression_scores0, w0)
    #sim_generations(population0, w0)

    print("Done gene{}.".format(genenumber))









#RUNNING THE PROGRAM
genenumber = 1
run_simulator(tfbs1, tfbs_len1, CDS_len1)
genenumber = 2
run_simulator(tfbs2, tfbs_len2, CDS_len2)