## Toolkit that takes in per site or windowed data and loads it in as a dataframe
## The class has a method to annotate SNPs and windows to genes if you load in a gene map
## The class has a method to annotate genes with Gene Ontology terms when gene lists are loaded
## The class has meethods to perform a Z test to compare a subset of the data to the entire dataset for specific GO terms



#importing modules
import os,sys
import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev
import scipy.stats
import matplotlib.pyplot as plt


#loci table class for anotating data such as Fst or pi and finding statistical enrichment for genes
class loci_table():

    def __init__(self, datafile, locitype="sites", win_size=0):
        """
        This class is currently configured to the Dmel genome from Flybase release 6.

        In order to use the .mapped2genes() method, this directory needs gene map csv files for each Chromosome
        In onrder to use the .annotate_GO() method, this directory must have GO gene lists in .txt format

        The loci_table class takes a csv or tsv file with the first column giving the chromosome and the second column giving the position.

        The class has two optional parameters. locitype denotes the kind of loci data, per site or windowed. locitype defaults to "sites" for windowed data pass "windowed"
        win_size denotes the step of the window. For per site data, this value defaults to 0. Pass an integer to change the window size to match your windowed data.
        """


        self.locitype = locitype
        self.win_size = win_size

        if isinstance(datafile,str) and ".csv" in datafile:
            self.table = pd.read_csv(datafile)
        elif isinstance(datafile,str) and ".pi" in datafile:
            self.table = pd.read_csv(datafile, sep='\t')
        else:
            print("ERROR: file not the correct type")

        self.ChrX_map = pd.read_csv('Dmel_genemap_ChrX.csv')
        self.Chr2L_map = pd.read_csv('Dmel_genemap_Chr2L.csv')
        self.Chr2R_map = pd.read_csv('Dmel_genemap_Chr2R.csv')
        self.Chr3L_map = pd.read_csv('Dmel_genemap_Chr3L.csv')
        self.Chr3R_map = pd.read_csv('Dmel_genemap_Chr3R.csv')

        #enrichment permutation dataframe
        self.permutation_df = pd.DataFrame()

        #top counts dictionary
        self.top_count = {}

        #string for summary printing
        self.summary = str(self.table.head(5)) + '\n.........................\n' + str(self.table.tail(5))

        print("####################################################################################################################################")
        print("This class is currently configured to the Dmel genome from Flybase release 6.")
        print("In order to use the .map2genes() method, this directory needs gene map csv files for each Chromosome")
        print("In order to use the .annotate_GO() method, this directory must have GO gene lists in .txt format")
        print("In order to use .enrichment_test(), plot_dist_param(), or .plot_enrichment_test() methods, the other methods must have been run")
        print("####################################################################################################################################\n\n")

        
    #CLASS FUNTIONS: You should not call these as they are used on the main methods of the class
    ##################################################################################################################
            
    #permutator methods
    def persite_null_permutator(self, n_permutes, nofset, term, GO_list):
        """
        This is the per site null permutator function.
        It randomly permutates through the table with a certain subset sample size and counts the number of genes represented by the loci and the weighted SNPs per genes (SNPs per kbp) for each permutation.
        It takes 3 parameters:
        n_permutes: integer, the number of random permutations of randomly subseting the dataframe
        nofset: integer, the size of the subset, get this value from the size of the top percentage of loci
        term: string, the GO term of interest
        """

        print("\n\npermutating...\n\n")
       
        #making lists
        gene_ls = self.ChrX_map['FBgn'].tolist() + self.Chr2L_map['FBgn'].tolist() + self.Chr2R_map['FBgn'].tolist() + self.Chr3L_map['FBgn'].tolist() + self.Chr3R_map['FBgn'].tolist()
        gene_len_ls = self.ChrX_map['length'].tolist() + self.Chr2L_map['length'].tolist() + self.Chr2R_map['length'].tolist() + self.Chr3L_map['length'].tolist() + self.Chr3R_map['length'].tolist()
        #making dictionaries
        genes_length = dict(zip(gene_ls, gene_len_ls))
        genes_length.update({'': 1})

        misc_df = pd.read_csv("misc_genes.csv")

        misc_gene_length = dict(zip(misc_df["FBgn"].to_list(), misc_df["Length"].to_list()))
        genes_length.update(misc_gene_length)




        #setting lists
        GO = []
        weighted_GO = []
        
        counter = 0
        #permutation
        for i in range(n_permutes):

            counter += 1
            if (counter % 100) == 0:
                print("Permutation:", str(counter), "of", n_permutes)
            
            #random sampling
            sampled_df = self.table.sample(n=nofset)


            #Neurogenesis
            GO_df = sampled_df.loc[sampled_df[term] == term]

            GO_gene_temporary = GO_df['FBgn'].tolist()
            GO_gene_temp = [i for i in GO_gene_temporary if i != "None"]
            GO_gene_str = "; ".join(GO_gene_temp)
            full_GO_gene_ls = GO_gene_str.split("; ")  ###

            GO_unique_gene_ls_temp = list(set(full_GO_gene_ls))
            GO_unique_gene_ls = [j for j in GO_unique_gene_ls_temp if j in GO_list]
            
            if '' in GO_unique_gene_ls:
                GO_unique_gene_ls.remove('') ###
            else:
                GO_unique_gene_ls = GO_unique_gene_ls

            #unweighted
            GO.append(len(GO_unique_gene_ls))

            #weighted
            GO_SNP_den = []

            for gene in GO_unique_gene_ls:
                gene_len = genes_length[gene]

                if len(GO_unique_gene_ls) == 0:
                    GO_SNP_den.append(0)
                else:
                    GO_SNP_den.append(((full_GO_gene_ls.count(gene)) / (gene_len)) * 1000)

            weighted_GO.append(sum(GO_SNP_den))
      


        #creating permutation dataframe
        permut_dict = {term:GO, term+'_weighted':weighted_GO}

        permut_df = pd.DataFrame.from_dict(permut_dict)

        return permut_df

    def win_null_permutator(self, n_permutes, nofset, term, GO_list):
        """
        This is the windowed null permutator function.
        It randomly permutates through the table with a certain subset sample size and counts the number of genes represented by the loci and the weighted SNPs per genes (SNPs per kbp) for each permutation.
        It takes 3 parameters:
        n_permutes: integer, the number of random permutations of randomly subseting the dataframe
        nofset: integer, the size of the subset, get this value from the size of the top percentage of loci
        term: string, the GO term of interest
        """

        print("\n\npermutating...\n\n")
       
        #making lists
        gene_ls = self.ChrX_map['FBgn'].tolist() + self.Chr2L_map['FBgn'].tolist() + self.Chr2R_map['FBgn'].tolist() + self.Chr3L_map['FBgn'].tolist() + self.Chr3R_map['FBgn'].tolist()
        gene_len_ls = self.ChrX_map['length'].tolist() + self.Chr2L_map['length'].tolist() + self.Chr2R_map['length'].tolist() + self.Chr3L_map['length'].tolist() + self.Chr3R_map['length'].tolist()
        



        #setting lists
        GO = []
        
        counter = 0
        #permutation
        for i in range(n_permutes):

            counter += 1
            if (counter % 1000) == 0:
                print("Permutation:", str(counter), "of", n_permutes)
            
            #random sampling
            sampled_df = self.table.sample(n=nofset)


            #Neurogenesis
            GO_df = sampled_df.loc[sampled_df[term] == term]

            GO_gene_temporary = GO_df['FBgn'].tolist()
            GO_gene_temp = [i for i in GO_gene_temporary if i != "None"]
            GO_gene_str = "; ".join(GO_gene_temp)
            full_GO_gene_ls = GO_gene_str.split("; ")
            GO_unique_gene_ls_temp = list(set(full_GO_gene_ls))
            GO_unique_gene_ls = [j for j in GO_unique_gene_ls_temp if j in GO_list]
            
            if '' in GO_unique_gene_ls:
                GO_unique_gene_ls.remove('') ###
            else:
                GO_unique_gene_ls = GO_unique_gene_ls

            #unweighted
            GO.append(len(GO_unique_gene_ls))

        #creating permutation dataframe
        permut_df = pd.DataFrame()

        permut_df[term] = GO

        return permut_df

    #making full and unique lists of genes from the table
    def top_gene_ls_maker(self, df):
        """
        makes a gene list and a unique gene list for the top percentage of data.
        Takes 1 parameter:
        df: dataframe of the top percentage
        Returns: full gene list and a unique gene list"""

        #making full gene list
        full_gene_temporary = df['FBgn'].tolist()
        full_gene_str = "; ".join(full_gene_temporary)
        full_gene_ls = full_gene_str.split("; ")    ####

        #making list of unique genes
        gene_ls = list(set(full_gene_ls))
        if "None" in gene_ls:
            gene_ls.remove("None") 
        else:
            gene_ls = gene_ls

        return full_gene_ls, gene_ls

    #Counting genes of a GO class
    def persite_top_GO_counter(self, full_ls, ls, GO_list):
        """
        This function counts the number of SNPs on the genes of interest and the SNP density of the genes of interest (SNP per kbp) of the top percentage.
        The function takes 3 parameters:
        full_ls: full list of genes with SNP hits in the top percentage
        ls: unique list of genes represented by the top percentage
        GO_list: list of genes of the GO class of interest
        """


        #making lists
        gene_ls = self.ChrX_map['FBgn'].tolist() + self.Chr2L_map['FBgn'].tolist() + self.Chr2R_map['FBgn'].tolist() + self.Chr3L_map['FBgn'].tolist() + self.Chr3R_map['FBgn'].tolist()
        gene_len_ls = self.ChrX_map['length'].tolist() + self.Chr2L_map['length'].tolist() + self.Chr2R_map['length'].tolist() + self.Chr3L_map['length'].tolist() + self.Chr3R_map['length'].tolist()
        #making dictionaries
        genes_length = dict(zip(gene_ls, gene_len_ls))
        genes_length.update({'': 1})

        misc_df = pd.read_csv("misc_genes.csv")

        misc_gene_length = dict(zip(misc_df["FBgn"].to_list(), misc_df["Length"].to_list()))
        genes_length.update(misc_gene_length)

        #counters for each GO term
        GO_count = 0
        weighted_GO = 0


        #neurogenesis loop
        for i in ls:
            if i in GO_list:
                GO_count += 1
                weighted_GO += (((full_ls.count(i)) / (genes_length[i])) * 1000)
            else:
                continue

            
        return GO_count, weighted_GO

    def win_top_GO_counter(self, ls, GO_list):
        """
        This function counts the number of genes of interest represented by the top percentage.
        The function takes 2 parameters:
        ls: unique list of genes represented by the top percentage
        GO_list: list of genes of the GO class of interest
        """

        #making lists
        gene_ls = self.ChrX_map['FBgn'].tolist() + self.Chr2L_map['FBgn'].tolist() + self.Chr2R_map['FBgn'].tolist() + self.Chr3L_map['FBgn'].tolist() + self.Chr3R_map['FBgn'].tolist()
        gene_len_ls = self.ChrX_map['length'].tolist() + self.Chr2L_map['length'].tolist() + self.Chr2R_map['length'].tolist() + self.Chr3L_map['length'].tolist() + self.Chr3R_map['length'].tolist()
        #making dictionaries
        genes_length = dict(zip(gene_ls, gene_len_ls))
        genes_length.update({'': 1})

        #counters for each GO term
        GO_count = 0


        #neurogenesis loop
        for i in ls:
            if i in GO_list:
                GO_count += 1
            else:
                continue

            
        return GO_count
    
    #Z-Test
    def ztest(self, topcount, nullcounts):
        """
        static function that performs a right tailed Z enrichment test between the set of interest and the null (whole) set
        The function takes 2 parameters:
        topcount: integer or float that represents the number of genes of interest in the set of interest
        nullcounts: list of integers or floats that represent the counts of the genes of interest in the null set from the permutation
        returns: z score and p value"""

        zscore = (topcount - mean(nullcounts)) / stdev(nullcounts)

        pvalue = scipy.stats.norm.sf(zscore)

        return zscore, pvalue

    ##################################################################################################################


    

    #gene mapper method
    def map2genes(self):
        """
        This method takes no arguments. It will map loci in the table to genes using the flybase genemap files.
        """

        print("mapping genes...")

        #getting list of coordinates
        coord_df = pd.DataFrame()
        coord_df['coord'] = self.table.iloc[:,0] + ':' + self.table.iloc[:,1].astype(str)
        
        self.coordinates = coord_df['coord'].to_list()

        

        #gene lists
        FBgn_list = []
        gene_symbol_list = []

        #looping through coordinates
        for i in self.coordinates:

            if 'X' in i:
                df =  self.ChrX_map
            elif '2L' in i:
                df =  self.Chr2L_map
            elif '2R' in i:
                df =  self.Chr2R_map
            elif '3L' in i:
                df =  self.Chr3L_map
            elif '3R' in i:
                df =  self.Chr3R_map
            else:
                print("Chromosome not found.")
                break
            
            #setting temp lists
            FBgn_list_temp = []
            gene_symbol_temp = []


            #looping for gene mapping (different for per site and windowed)
            if self.locitype == "windowed":

                for row in df.itertuples():
                    #tuple to list
                    row_list = list(row)

                    if ((int(i.split(':')[1]) >= row_list[4]) and (int(i.split(':')[1]) <= row_list[5])) and (((int(i.split(':')[1])+self.win_size)) >= row_list[4]) and (((int(i.split(':')[1])+self.win_size)) >= row_list[5]):
                        FBgn_list_temp.append(row_list[1])
                        gene_symbol_temp.append(str(row_list[2]))
                    elif ((int(i.split(':')[1]) <= row_list[4]) and (int(i.split(':')[1]) <= row_list[5])) and (((int(i.split(':')[1])+self.win_size)) >= row_list[4]) and (((int(i.split(':')[1])+self.win_size)) >= row_list[5]):
                        FBgn_list_temp.append(row_list[1])
                        gene_symbol_temp.append(str(row_list[2]))
                    elif ((int(i.split(':')[1]) <= row_list[4]) and (int(i.split(':')[1]) <= row_list[5])) and (((int(i.split(':')[1])+self.win_size)) > row_list[4]) and (((int(i.split(':')[1])+self.win_size)) < row_list[5]):
                        FBgn_list_temp.append(row_list[1])
                        gene_symbol_temp.append(str(row_list[2]))
                    else:
                        continue
                
                #temp lists to strings
                FBgn_str = "; ".join(FBgn_list_temp)
                gene_str = "; ".join(gene_symbol_temp)

                #checking if there is no hit
                if len(FBgn_str) == 0:
                    FBgn_list.append('None')
                    gene_symbol_list.append('None')
                else:
                    FBgn_list.append(FBgn_str)
                    gene_symbol_list.append(gene_str)


            
            else:

                for row in df.itertuples():
                    #tuple to list
                    row_list = list(row)

                    if int(i.split(':')[1]) >= row_list[4] and int(i.split(':')[1]) <= row_list[5]:
                        FBgn_list_temp.append(row_list[1])
                        gene_symbol_temp.append(str(row_list[2]))
                    else:
                        continue
                
                #temp lists to strings
                FBgn_str = "; ".join(FBgn_list_temp)
                gene_str = "; ".join(gene_symbol_temp)

                #checking if there is no hit
                if len(FBgn_str) == 0:
                    FBgn_list.append('None')
                    gene_symbol_list.append('None')
                else:
                    FBgn_list.append(FBgn_str)
                    gene_symbol_list.append(gene_str)



        #making dataframe
        dictionary = {'FBgn' : FBgn_list, 'gene_symbol' : gene_symbol_list}

        gene_hits_df = pd.DataFrame(dictionary)

        #merging with self.table
        self.table = pd.merge(self.table, gene_hits_df, left_index=True, right_index=True)

        #string for summary printing
        self.summary = str(self.table.head(5)) + '\n.........................\n' + str(self.table.tail(5)) + '\n\n########################################################################################################################################################################################################################'




    #GO annotator method
    def annotate_GO(self, GO_listfile, GO_term):
        """
        This method takes two parameters:
        GO_listfile: string of a filename that contains a list of genes, one gene per line, that fir your chosen GO term
        GO_term: string of the GO term you whan to annotate.
        """

        print("annotating genes with GO term")

        ##Flybase GO lists
        FB_GO_file = open(GO_listfile, "r")
        FB_GO_genes = FB_GO_file.read()
        #to list
        self.FB_GO_list = FB_GO_genes.split("\n")
        FB_GO_file.close()

        
        #gene list from self.table
        gene_list = self.table['FBgn'].tolist()

       

        #setting list
        GO = []

        #looping through null genes
        for gene in gene_list:
            gene_ls = gene.split("; ")

            str_GO = "None"
            
            for i in gene_ls:
                if i in self.FB_GO_list:
                    str_GO = GO_term
                else:
                    continue


            GO.append(str_GO)
            

        #annotating self.table
        dictionary = {GO_term : GO}
        GO_df = pd.DataFrame(dictionary)

        self.table = pd.merge(self.table, GO_df, left_index=True, right_index=True)

        #string for summary printing
        self.summary = str(self.table.head(5)) + '\n.........................\n' + str(self.table.tail(5)) + '\n\n########################################################################################################################################################################################################################'

        
 

    def enrichment_test(self, numofpermutes, parameter, top_percentage, GO_term, GO_listfile):
        """
        This method performs a enrichment test (z-test) between the top percentage of loci (based on your data's parameter) and the null (whole) set.
        This method takes 3 parameters:
        paramter: string of column name in the table that you want to take the top percentage of (ex: 'Avg_Fst')
        top_percentage: an integer or float in units of percent that denotes what percentage of the table's data should be used as the set of interest
        GO_term: string of the GO term of interest
        """

        ##Flybase GO lists
        FB_GO_file = open(GO_listfile, "r")
        FB_GO_genes = FB_GO_file.read()
        #to list
        FB_GO_ls = FB_GO_genes.split("\n")
        FB_GO_file.close()

        #size of top percentage df
        top_len = int(len(self.table) * (top_percentage/100))

        if top_len < 10:
            print('\n############################################################################################################')
            print("WARNING: Percentage threhold only represents", top_len, "loci. Consider using higher percentage threshold.")
            print('\n############################################################################################################')

        #top percentage df
        top_df = self.table.sort_values(by=parameter, ascending=False)
        top_df = top_df.iloc[0:top_len]


        if self.locitype == "windowed":

            #top counting
            full, genes = self.top_gene_ls_maker(top_df)
            top_count = self.win_top_GO_counter(genes, FB_GO_ls)

            self.top_count.update({GO_term: top_count})

            #null permutation
            perm_df = self.win_null_permutator(numofpermutes, top_len, GO_term, FB_GO_ls)

            if self.permutation_df.empty:
                self.permutation_df = perm_df
            else:
                self.permutation_df = pd.merge(self.permutation_df, perm_df, left_index=True, right_index=True)


            null_counts = perm_df.iloc[:, 0]

            #z-test adn p-value
            z, p = self.ztest(top_count, null_counts)

            string = '\n\n\n\n' + GO_term + ' Top ' + str(top_percentage) + '% ' + 'of loci based on highest ' + str(parameter) + '\n' + 'z: ' + str(z) + ' p: ' + str(p) + '\n\n'
            self.summary += string
            print(string)


        else:

            #top counting
            full, genes = self.top_gene_ls_maker(top_df)
            top_count, top_weighted_count = self.persite_top_GO_counter(full, genes, FB_GO_ls)

            self.top_count.update({GO_term: top_count, GO_term+'_weighted': top_weighted_count})

            #null permutation
            self.permutation_df = self.persite_null_permutator(numofpermutes, top_len, GO_term, FB_GO_ls)

            null_counts = self.permutation_df.iloc[:, 0]
            null_weighted_counts = self.permutation_df.iloc[:, 1]

            #z-test adn p-value
            z, p = self.ztest(top_count, null_counts)

            zweight, pweight = self.ztest(top_weighted_count, null_weighted_counts)

            string = '\n\n\n\n' + GO_term + ' Top ' + str(top_percentage) + '% ' + 'of loci based on highest ' + str(parameter) + '\n' + 'z: ' + str(z) + ' p: ' + str(p) + '\n' + 'z weighted: ' + str(zweight) + ' p weighted: ' + str(pweight) + '\n\n'
            self.summary += string
            print(string)

            

    def plot_dist_param(self, parameter):
        """
        This function plots the distribution of the paramer of your choosing in the table.
        This function takes 1 parameter:
        parameter: string of column name over which to plot the distribution.
        """

        param_ls = self.table[parameter].to_list()
        param_a = np.array(param_ls)
       
        
        #number of bins, Freedman–Diaconis
        q25, q75 = np.percentile(param_a, [25, 75])
        bin_width = 2 * (q75 - q25) * len(param_a) ** (-1/3)
        bins = round((param_a.max() - param_a.min()) / bin_width)
        
        
        
        #plotting
        plt.hist(param_a, density=True, bins=bins )
        plt.ylabel('Density')
        plt.xlabel(parameter)
        plt.title(parameter + " Distribution")

        mn, mx = plt.xlim()
        plt.xlim(mn, mx)
        kde_xs = np.linspace(mn, mx, 300)
        kde = scipy.stats.gaussian_kde(param_a)
        plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")

        plt.savefig(parameter + '.png', bbox_inches='tight')

        plt.clf()  




    def plot_enrichment_test(self, plotname, term, weighted=False):
        """,
        This function is for plotting the distribution of the null counts compared to the count of the top percentage from the last run enrichment test.
        This function takes 3 parameters:
        plotname: str, name of file to save the figure as
        plottitle: str, title of the plot
        term: string of the GO term of interest
        weighted: bool, optional. whether or not to select the counts permutation or the weighted permutation for per site data. Defaults to False
        """

        if weighted == True:
            colname = term + '_weighted'
            word = ' SNP density (SNPs per kbp)'
        else:
            colname = term
            word = ' Gene Counts'

        top = self.top_count[colname]
        null = np.array(self.permutation_df[colname].to_list())

        
        #number of bins, Freedman–Diaconis
        q25, q75 = np.percentile(null, [25, 75])
        bin_width = 2 * (q75 - q25) * len(null) ** (-1/3)

        if bin_width == 0:
            bins = 10
        else:
            bins = round((null.max() - null.min()) / bin_width) 

        
        
        #plotting
        plt.hist(null, density=True, bins=bins)
        plt.axvline(x=top, color='r', lw=3)
        plt.ylabel('Density')
        plt.xlabel(term + word)
        #plt.title(plottitle)

        mn, mx = plt.xlim()
        plt.xlim(mn, mx)
        kde_xs = np.linspace(mn, mx, 300)
        kde = scipy.stats.gaussian_kde(null)
        plt.plot(kde_xs, kde.pdf(kde_xs), label="PDF")


        plt.savefig(plotname + '.png', bbox_inches='tight')

        plt.clf()  



        
    def __str__(self):
        """Overload of str for printing and saving to file"""

        return self.summary + '\n\nTop Counts:\n' + str(self.top_count) + '\n\n\n'

            

        

















    













    


