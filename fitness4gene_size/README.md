## Gene Expression, Fitness, and Evolution Simulator

This program simulates the evolution of gene expression and fitness in a Wright-Fisher model using a diploid population of obligate outbreeding hermaphrodites. The goal of the simulation is to compare two differently sized genes, with one gene having a longer coding sequence (CDS) and more transcription factor binding sites (TFBS) than the other. The simulation parameters are set through command-line arguments.

### Installation

Clone this repository, navigate to the directory, and create a conda environment for the tool as follows:

```
conda env create -f FITNESS_SIM.yml
```

### Usage

Run the script through a terminal using a command such as the following:

```
python fitness4gene_size.py --mu 0.001 --pop 100 --gen 100 --express2fit linear --burnin 50 --seed 1234 --tfbs1 10 --tfbslen1 20 --cds1 2000 --tfbs2 20 --tfbslen2 30 --cds2 3000
```

### Required Arguments

The script requires the following arguments:

- `--express2fit`: A string specifying the fitness function to be used. The options are "linear", "parabolic", or "sigmoidal". Example: `--express2fit linear`.
- `--seed`: An integer specifying the seed value for randomization modules. Example: `--seed 1234`.
- `--tfbs1`, `--tfbslen1`, and `--cds1`: Integers specifying the number of TFBSs, the average length of TFBSs, and the length of the CDS for gene 1, respectively. Example: `--tfbs1 10 --tfbslen1 20 --cds1 2000`.
- `--tfbs2`, `--tfbslen2`, and `--cds2`: Integers specifying the number of TFBSs, the average length of TFBSs, and the length of the CDS for gene 2, respectively. Example: `--tfbs2 20 --tfbslen2 30 --cds2 3000`.

Note that starting in version 1.7, the tfbs, tfbslen, and cds arguments are now technically optional to allow for running the simulation for only one gene (either gene 1 or gene 2).  
However, you must set either tfbs1, tfbslen1, and cds1; tfbs2, tfbslen2, and cds2; or all of them.  

### Optional Arguments

The script also accepts the following optional arguments, which have default values:

- `--mu`: A float specifying the per-site mutation rate. The default value is 0.001. Example: `--mu 0.002`.
- `--pop`: An integer specifying the constant population size. The default value is 100. Example: `--pop 50`.
- `--gen`: An integer specifying the number of generations for the simulation. The default value is 100. Example: `--gen 200`.
- `--burnin`: An integer specifying the number of burn-in generations. The default value is 50. Example: `--burnin 100`.
- `--mueffects`: A string specifying the effect of mutations to be used as a control. The default is "all", change to "deleterious" for all mutations to be deleterious. Example: `--mueffects deleterious`
- `--mutypes`: A string specifying if mutations land in TFBSs, CDS, or both to be used as a control. The default is "all", The other options are "CDS" or "CRM". Example: `--mutypes CRM`
- `--outputdir`: A string specifying the path to the directory in which to save the output files, note: the directory sting must end in "/" and it can be a relative path if the output directory is a subdirectory to the one in which this script is saved.. Example: `--outputdir path/to/save/output/`



## Model

In this simplified model, gene expression occurs when a transcription factor (TF) binds to at least one TF binding site (TFBS) within the cis-regulatory module (CRM). Longer coding sequences (CDS) have a higher likelihood of experiencing non-synonymous mutations.

Each gene contains a specific number of TFBS, each with a probability score indicating its likelihood of being bound by a TF. The cis-regulatory model score represents the probability of at least one TFBS being bound at any given time. The CDS initially receives a score of 1, indicating a normal gene product, with lower scores indicating reduced functionality and higher scores indicating increased functionality.

During each generation, mutations are introduced based on a per-site mutation rate, sampled from a binomial distribution. Mutations occurring within TFBS affect the TFBS score proportionally. Mutation effects are sampled from a gamma distribution, where approximately 10% of mutations increase binding probability, while on average, the binding affinity decreases to around 60% of the original affinity.

Mutations within the CDS affect gene function only if they are non-synonymous, which accounts for approximately 73% of all CDS mutations. These mutations also impact function proportionally. Mutational effect proportions are sampled from a gamma distribution, where non-synonymous mutations are expected to decrease gene expression by approximately 65%, with around 1% of mutations having at least a neutral effect.

Cis-regulatory model scores reflect relative expression levels and are input into one of three fitness functions for selection at each generation: linear, parabolic, or sigmoidal. The CRM's fitness contribution score is computed based on the chosen fitness function. The gene score is then determined by multiplying the CRM fitness score and the CDS score. The gene score is calculated by multiplying the cis-regulatory model score and the CDS score. As the species is diploid, the scores of both homologs are averaged to obtain an individual score.

The population consists of diploid individuals that are obligate outbreeding hermaphrodites, and the generation size remains constant.