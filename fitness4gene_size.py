"""
fitness4gene_size version 1.2

Author: Phil Baldassari

Description: Wright-Fisher process on two differently sized genes. One gene has a longer CDS and more TFBS than the other.
In this simplistic model, a gene is expressed when each TFBS is bound by a TF at least once during a specified time period in ms.
A longer CDS has more chance for Non-synonymous mutation.

Model:
Each gene has a certain number of TFBS each with a probability score of it being bound by a TF at least once during a specified time period in ms.
The Cis-regulatory model score cooresponds to the probability that each TFBS was bound at least once during the specified time period.
The CDS recieves an initial score of 1 meaning a normal gene product.
A score of 1 indicates a normal gene product with lower scores being less expressed and higher scores being more expressed.
Each generation mutations are added based on a per site mutation rate. the number of mutations are sampled from a binomial distribution.
Mutations in TFBS affect the TFBS score propotianally. Multiplpiers are smalled from a normal distribution with a mean of -0.01 which indicates 
that on average mutation decreases the expression contribution of that TFBS by 1%.
Mutations on CDS only affect function if they are nonsynonymous which occurs with around 73% of all mutations in CDS. These muations also affect function propotianlly.
Mutational effect proportions are sampled from a gamma distribution of shape 3 and rate 8.4 (scale 0.11905) so that the expected (average) effect of a NS mutation decreases the gene expression
by ~64% and mutations will be at least nearly neutral ~2% of the time and at least truly neutral ~1% of the time.
Effect propotional of multiple NS mutations are mutliplied together before being multiplied to the original expression level score.
The gene score is computed by multiplying Cis-regulatory model score and the CDS score. Since species is diploid, the scores of both homologs are averaged for an individual score.
These scores are analogous to a relative expression level. The scores are input into one of 3 fitness functions for selection at each genreation: i) linear, ii) parabolic, iii) sigmoidal.
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
parser.add_argument('--bindtime', type=float, default=10, help='Timeframe for TF binding in milliseconds, default is 10')
parser.add_argument('--express2fit', type=str, required=True, help='function repressenting relationship b/n expression and fitness, options are: "linear", "parabolic", or "sigmoidal"')
parser.add_argument('--burnin', type=int, default=50, help='number of burnin generations, default is 50')
parser.add_argument('--burninmode', type=str, default="regenerative", help='burnin mode, deafault is "regenerative" alternative mode is "neutral"')
parser.add_argument('--seed', type=int, required=True, help='seed value for randomization modules')

parser.add_argument('--tfbs1', type=int, required=True, help='number of TFBSs for gene 1')
parser.add_argument('--tfbslen1', type=int, required=True, help='average length of TFBSs for gene 1')
parser.add_argument('--cds1', type=int, required=True, help='length of CDS for gene 1 note: this number does NOT have to be divisible by 3')

parser.add_argument('--tfbs2', type=int, required=True, help='number of TFBSs for gene 2')
parser.add_argument('--tfbslen2', type=int, required=True, help='average length of TFBSs for gene 2')
parser.add_argument('--cds2', type=int, required=True, help='length of CDS for gene 2 note: this number does NOT have to be divisible by 3')

#parse args
args = parser.parse_args()

#assign variables
mu = args.mu
n = args.pop
g = args.gen
bindtime = args.bindtime
fitness_function = args.express2fit
burnin_g = args.burnin
burnin_mode = args.burninmode
seed = args.seed

tfbs1 = args.tfbs1
tfbs_len1 = args.tfbslen1
CDS_len1 = args.cds1

tfbs2 = args.tfbs2
tfbs_len2 = args.tfbslen2
CDS_len2 = args.cds2


#fitness function parameter
a = 1

#fitness functions
def linear(x):
    """f(x) = ax"""
    y = a * x

    return y

def parabolic(x):
    """f(x) = -(2x + a)^2 + 1"""

    y = -1 * ((2*x)+a)**2 + 1

    return y

def sigmoidal(x):
    """f(x) = 1/(1 + e^-5(2x - a))"""

    y = 1 / (1 +(math.e ** (-5*(2*x + a))))

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
    Function takes the list of TFBS scores and uses them to compute the probability of each TFBS being bound at least once in the set time period.
    Pr(e) = 1 - PI(1-b)^t
    Returns the numerical value
    """
    score = 1

    for i in ls:
        score *= (1 - i) ** bindtime

    final_score = 1 - score

    return final_score


def compute_indv_scores(pop_ls):
    """
    Function takes in a population list and computes a raw score for each individual by coputing a score for each gene and averaging the score for each gene pair (diploid individual)
    Scores are coputed by the product of all TFBS scores and the CDS score. The scores for each gene in the homologous pair are averaged adn the average is appended to the raw score list.
    Returns a list of raw scores for every indidivual in the population. This list must be normalized later in order to compute fitness.
    """

    indv_scores = []
    
    #looping through individuals
    for indv_idx in range(len(pop_ls)):
        gene_scores = []

        #loping through each gene per individual
        for gene_idx in range(2):
            #computing score for each gene
            score = compute_CRM_scores(pop_ls[indv_idx][gene_idx][0][0]) * pop_ls[indv_idx][gene_idx][1][0]
            gene_scores.append(score)

        #averaging score for each diploid pair
        indv_scores.append(mean(gene_scores))

    return indv_scores



#simulation functions
def burnin(num_of_tfbs, length_of_tfbs, length_of_cds):
    """
    Function to run a burnin of a set amount of generations (does not need to exceed 200). The expression scores are compiled in a large list.
    The average of the expression scores will be used to calibrate the fitness functions to set the average to a fitness of 0.5
    The final population is then returned to be used in the simulation used as the starting population.
    The burnin can run in two modes: regenerative or neutral. Regenerative will randomly generate a new starting population for each generation while
    the neutral mode runs a neutral W-F process. In both cases, scores are appended to a list for averaging.
    Returns starting population and the fitness function parameter.
    In order for this function to work properly fitnesses must be computed outside the function using the modifided fitness functions.
    """

    #setting parameters
    tfbs = num_of_tfbs
    tfbs_len = length_of_tfbs
    CDS_len = length_of_cds

    #burnin scores
    scores_ls = []

    #running burnin in either neutral or regenerative mode
    if burnin_mode == "neutral":
        populationBI0, scBI0, wBI0 = starting_pop(tfbs, tfbs_len, CDS_len)
        pplnBI = populationBI0
        wBI = wBI0

        scores_ls += scBI0

        for gen in range(burnin_g):
            pplnBI, scBI, wBI = next_gen(pplnBI, wBI, selection="neutral")

            scores_ls += scBI

    else:
        for gen in range(burnin_g):
            pplnBI, scBI, wBI = starting_pop(tfbs, tfbs_len, CDS_len)

            scores_ls += scBI

    #parameterizing fitness functions
    s_bar = mean(scores_ls)

    if fitness_function == "linear":
        param = 1/(2*s_bar)
    elif fitness_function == "parabolic":
        paramls = []
        param1 = (math.sqrt(0.5)) - (2*s_bar)
        param2 = -1 * ((math.sqrt(0.5)) + (2*s_bar))
        paramls.append(param1)
        paramls.append(param2)
        param = min(paramls)
    elif fitness_function == "sigmoidal":
        param = -2 * s_bar
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
    Returns: population list and list of scores
    """

    popls = []

    #tfbs lengths
    tfbs_lens = list(np.random.gamma((tfbs_len/2), scale=2, size=tfbs))
    tfbs_lens = [int(x) for x in tfbs_lens]

    #genreating genes for each individual
    for i in range(n):
        indv = []
        #each individual is diploid
        for j in range(2):
            gene = []

            #scores for TFBS binding probabilites (Lahdesmaki et al. 2008)??
            #tfbs_probs = list(np.random.exponential(0.1, size=tfbs))
            #getting rid of initial probabilities at zero
            #tfbs_probs = [p if p > 0 else 0.0001 for p in tfbs_probs]

            tfbs_probs = [0.1 for i in range(tfbs)]


            """
            #however, empirical data shows a second mode at 1 with frequency of 0.01 (Lahdesmaki et al. 2008)
            numof1s = np.random.binomial(tfbs, 0.01)
            where1s = random.choices([idx for idx in range(tfbs)], k=numof1s)

            for idx in where1s:
                tfbs_probs[idx] = 1
            """


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
    scores = compute_indv_scores(popls)

    #fitness effects
    fitnesses = []
    for s in scores:
        fitnesses.append(fitness_func[fitness_function](s))

    #bounding at zero
    scores = bound_at_zero(scores)
    fitnesses = bound_at_zero(fitnesses)

    return popls, scores, fitnesses


def next_gen(pop_ls, ws, selection="selection"):
    """
    Fuction outputs the next generation from a previous genreation using the W-F model.
    Inputs are the population list and the list of fitness scores.
    The fitness will be used to simulate mating with a skew toward more mating instances between individuals of higher fitness (selection).
    Returns a new population list as well as a new fitness score list.
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


        #adding mutations in TFBS which change the scores based on a normally distributed differential
        for gene in offspring_indv:
            #looping through TFBS
            for idx in range(len(gene[0][0])):
                #how many mutations for each TFBS
                howmanyTFBS = np.random.binomial(gene[0][1][idx], mu)

                #changing scores proportinally: s' = s(1+x)
                #pc_diff = 0
                pc_prop = 1
                for j in range(howmanyTFBS):
                    #pc_diff += np.random.normal(-0.01, 0.02)
                    pc_prop += np.random.gamma(3, 0.133333)
                
                if howmanyTFBS == 0:
                    continue
                else:
                    #new score making sure to keep within the bounds [0,inf)
                    #newscore = gene[0][0][idx] * (1 + pc_diff)
                    newscore = gene[0][0][idx] * pc_prop
                    """
                    if newscore < 0:
                        newscore = 0
                    else:
                        newscore = newscore
                    """

                    #setting new score
                    gene[0][0][idx] = newscore
            
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

            #changing scores proportinally: s' = s(1+x)
            cds_pc_prop = 1
            for j in range(SNP_eff.count("NS")):
                #cds_pc_prop *= np.random.gamma(3, 0.11905)
                cds_pc_prop *= np.random.gamma(2, 0.15)
            
            if howmanyCDS == 0:
                continue
            else:
                #new score
                cds_score = gene[1][0] * cds_pc_prop

            gene[1][0] = cds_score

        
        #adding offspring to population
        new_pop.append(offspring_indv)

    #computing raw scores for each individual
    scores = compute_indv_scores(new_pop)

    #fitness effects
    fitnesses = []
    for s in scores:
        fitnesses.append(fitness_func[fitness_function](s))

    #bounding at zero
    scores = bound_at_zero(scores)
    fitnesses = bound_at_zero(fitnesses)

    return new_pop, scores, fitnesses


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
    max_Dw = [0]
    avg_Dw = [0]
    var_Dw = [0]
    max_s.append(max(scores0))
    avg_s.append(mean(scores0))
    var_s.append(variance(scores0))
    avg_w.append(mean(fitnesses0))
    max_w.append(max(fitnesses0))
    var_w.append(variance(fitnesses0))



    #holding list for relative fitnesses for computing fitness deltas between generations
    w_previous = fitnesses0

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

        #computing fitness deltas
        Dw = np.array(w) - np.array(w_previous)
        max_Dw.append(Dw.max())
        avg_Dw.append(Dw.mean())
        var_Dw.append(Dw.var())

        #setting previous generation fitnesses
        w_previous = w




    #plotting
    fig, axs = plt.subplots(6, 1, figsize=(10, 15))

    axs[0].plot(generation, max_s)
    axs[0].set_title("Maximum Expression Score per Generation")
    axs[0].set_xlabel("Generation")
    axs[0].set_ylabel("Max Expression Score")

    axs[1].plot(generation, avg_s)
    axs[1].set_title("Average Expression Score per Generation")
    axs[1].set_xlabel("Generation")
    axs[1].set_ylabel("Avg Expression Score")

    axs[2].plot(generation, max_w)
    axs[2].set_title("Maximum Fitness per Generation")
    axs[2].set_xlabel("Generation")
    axs[2].set_ylabel("Max Relative Fitness")

    axs[3].plot(generation, avg_w)
    axs[3].set_title("Average Fitness per Generation")
    axs[3].set_xlabel("Generation")
    axs[3].set_ylabel("Avg Relative Fitness")

    axs[4].plot(generation, max_Dw)
    axs[4].set_title("Maximum Delta Fitness per Generation")
    axs[4].set_xlabel("Generation")
    axs[4].set_ylabel("Max Delta Fitness")

    axs[5].plot(generation, avg_Dw)
    axs[5].set_title("Average Delta Fitness per Generation")
    axs[5].set_xlabel("Generation")
    axs[5].set_ylabel("Avg Delta Fitness")

    fig.subplots_adjust(hspace=0.5)

    plt.savefig('WF_plot_gene{}.png'.format(genenumber))

    #saving csv
    dictionary = {
        "generation":generation, 
        "max_expression_score": max_s, "avg_expression_score": avg_s, "variance_expression_score":var_s, 
        "max_fitness": max_w, "avg_fitness": avg_w, "variance_fitness":var_w, 
        "max_DELTAfitness": max_Dw, "avg_DELTAfitness": avg_Dw, "variance_DELTAfitness":var_Dw
        }
    df = pd.DataFrame(dictionary)
    df.to_csv("WF_data_gene{}.csv".format(genenumber), index=False)


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

    #computing raw scores for each individual
    scores0 = compute_indv_scores(population0)
    #fitness effects
    fitnesses0 = []
    for s in scores0:
        fitnesses0.append(fitness_func[fitness_function](s))
    #bounding at zero
    sc0 = bound_at_zero(scores0)
    w0 = bound_at_zero(fitnesses0)

    print("Simulating generations...")
    #simulating generations
    sim_generations(population0, sc0, w0)

    print("Done gene{}.".format(genenumber))




#RUNNING THE PROGRAM
genenumber = 1
run_simulator(tfbs1, tfbs_len1, CDS_len1)
genenumber = 2
run_simulator(tfbs2, tfbs_len2, CDS_len2)