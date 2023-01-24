#Cluster-gene association comparator
import os
import re
from collections import defaultdict

def genelistreader(genelistdir):
    genome_to_namedgene = defaultdict(list)
    gene_to_IGRtagGENOME = defaultdict(list)
    genelist_list = os.listdir(genelistdir) #I am a master of naming variables
    print(f"Beginning Genelist readings ({len(genelist_list)} total)")
    for genelist in genelist_list:
        beginreading = False
        id = re.findall("GENELIST_(.*)\.txt",genelist)[0]
        with open(str(genelistdir+"\\"+genelist), "r") as file:
            for line in file:
                if beginreading == True:
                    splitline = line.split()
                    if len(splitline) == 4:
                        gene_to_IGRtagGENOME[splitline[2]].append([splitline[3],splitline[1],id]) #key = gene name, value = [clusterID,locustag,genomeID]
                        if id in genome_to_namedgene.keys():
                            if splitline[2] not in genome_to_namedgene[id] and splitline[2] != "Unknown_Gene":
                                genome_to_namedgene[id].append(splitline[2])
                        elif id not in genome_to_namedgene.keys():
                            if splitline[2] != "Unknown_Gene":
                                genome_to_namedgene[id].append(splitline[2])
                    else:
                        print(f"Potential error in {genelist} at line {line}")
                if "LOCUS TAG TO GENE LIST" in line:
                    beginreading = True
    print("Genelists read, moving to association comparisons")
    return gene_to_IGRtagGENOME, genome_to_namedgene

def association_comparator(genetoIGG): #IGG shorthand for IGR gene genome, variable is the dict produced by the above function
    significant_gene_freq_dict = defaultdict(list)
    gene_freq_dict = defaultdict(list)
    for gene in genetoIGG.keys():
        if gene != "Unknown_Gene":
            #print(f"Beginning analysis of {gene}")
            singleclusteroccs = [] #for each gene, holds all clusterIDs associated with it WITHOUT DUPLICATES - Used to tell .count() function what to count, without entering same clusterID twice
            totalclusteroccs = [] #for each gene, holds all clusterIDs associated with it WITH DUPLICATES - For clusterID counting purposes
            cluster_count_list = [] #nested list, each entry in format [cluster, numberoccurrencesofcluster,percentoftotaloccs]
            for cluster,tag,genome in genetoIGG[gene]:
                totalclusteroccs.append(cluster)
                if cluster not in singleclusteroccs:
                    singleclusteroccs.append(cluster)
            denominator = len(totalclusteroccs)
            expected_if_random_dist = ((denominator/len(singleclusteroccs))/denominator)*100
            for cluster in singleclusteroccs:
                numberassocs = totalclusteroccs.count(cluster)
                fractionassoc = (numberassocs/denominator)*100
                cluster_count_list.append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist]) 
                if fractionassoc > (15+expected_if_random_dist) and len(singleclusteroccs) != 1: #15 is arbitrary, checks if fraction assoc is significantly higher than expected.
                    print(f"{cluster} found associated with {numberassocs} copies of {gene} which occurred {denominator} times across all genomes")
                    significant_gene_freq_dict[gene].append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist]) 
                if len(singleclusteroccs) == 1:
                    print(f"{cluster} occurs {numberassocs} times and IS THE ONLY IGR ASSOCIATED WITH {gene} which occurred {denominator} times across all genomes")
            #print("############################################################################################")
            gene_freq_dict[gene].extend(cluster_count_list)
    return significant_gene_freq_dict, gene_freq_dict



validchoice = False
while validchoice == False:
    datasetprompt = input("What dataset are you wanting to list genes for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        genelistloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_both_genomes"
        outputloc = "Thesis_project_data\\Data\\Processed_data\\IGR_association_analysis_both_genomes"
        mastloc = "Thesis_project_data\\Data\\Processed_data\\mastout"
        mergedclustersloc = "Thesis_project_data\\Data\\Processed_data\\mast_merged_igrs_bothgenomes"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        genelistloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_comp_genomes"
        outputloc ="Thesis_project_data\\Data\\Processed_data\\IGR_association_analysis_comp_genomes"
        mastloc = "Thesis_project_data\\Data\\Processed_data\\mastout"
        mergedclustersloc = "Thesis_project_data\\Data\\Processed_data\\mast_merged_igrs_bothgenomes"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")

unformattedids = os.listdir(genelistloc)

gene_to_igrgenegenome_dict,genome_to_gene = genelistreader(genelistloc)

sig_gene_dict, gene_dict = association_comparator(gene_to_igrgenegenome_dict)
print("Writing the collated file summarising all IGR-gene associations across all genomes")
with open(str(outputloc+"\\Collated_gene_igr_associations.txt"), "w") as file:
    file.write("SIGNIFICANT GENE CLUSTER ASSOCIATIONS OCCURRING >100x WITH FRACTION ASSOCIATED +15% HIGHER THAN EXPECTED\n")
    file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
    for gene in sig_gene_dict.keys():
        for cluster,denominator,numberassoc,fracassoc,expected_dist in sig_gene_dict[gene]:
            if numberassoc > 100:
                file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{fracassoc}\t{expected_dist}\n")
                print(f"Look into {gene} downstream of {cluster} a total of {numberassoc} times ({fracassoc} of all occurrences, expected {expected_dist})")
    file.write("\n#######################################################################\n")
    file.write("ALL GENE TO IGR ASSOCIATIONS ACROSS ALL GENOMES\n")
    file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
    for gene in gene_dict:
        for cluster,denominator,numberassoc,fracassoc,expected_dist in gene_dict[gene]:
            file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{fracassoc}\t{expected_dist}\n")


for ufid in unformattedids:
    id = re.findall("GENELIST_(.*)\.txt",ufid)[0]
    print(f"Writing gene-igr associations for {id}")
    with open(str(outputloc+"\\IGR_associations_pergenome\\Assoc_analysis_"+id+".txt"),"w") as file:
        if id in genome_to_gene.keys():
            file.write(f"{id} GENE TO IGR ASSOCIATIONS\n")
            file.write("SIGNIFICANT GENE CLUSTER ASSOCIATIONS OCCURRING >100x WITH FRACTION ASSOCIATED +15% HIGHER THAN EXPECTED\n")
            file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
            for gene in genome_to_gene[id]:
                clusteroccs = []
                for igr,tag,genome in gene_to_igrgenegenome_dict[gene]:
                    if genome == id:
                        clusteroccs.append(igr)
                for cluster,denominator,numberassoc,fracassoc,expected_dist in sig_gene_dict[gene]:
                    if cluster in clusteroccs:
                        file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{fracassoc}\t{expected_dist}\n")
            file.write("\n#######################################################################\n")
            file.write(f"ALL GENE TO IGR ASSOCIATIONS ACROSS {id}\n")
            file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
            for gene in genome_to_gene[id]:
                clusteroccs = []
                for igr,tag,genome in gene_to_igrgenegenome_dict[gene]:
                    if genome == id:
                        clusteroccs.append(igr)
                for cluster,denominator,numberassoc,fracassoc,expected_dist in gene_dict[gene]:
                    if cluster in clusteroccs:
                        file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{fracassoc}\t{expected_dist}\n")
        elif id not in genome_to_gene.keys():
            print(f"An error has occurred, {id} doesn't have any listed named genes")