#MAST_contents_page_maker
import os
from collections import defaultdict
import re
import networkx as nx


def mastreader(mastloc,genomeids):
    genomes_cluster_dict = defaultdict(list) #For each cluster, returns the genome ID(s) of the genomes for which the mast.txt contained the cluster
    scoringseqfound = False
    for id in genomeids:
        print(f"Beginning reading of MAST.txt for {id}")
        masttxt = str(mastloc+"\\"+id+"\\mast.txt")
        with open(masttxt, "r") as file:
            for line in file:
                if "SECTION I: HIGH-SCORING SEQUENCES" in line:
                    scoringseqfound = True
                if "SECTION II: MOTIF DIAGRAMS" in line:
                    scoringseqfound = False
                    break
                if scoringseqfound == True:
                    if "Cluster_" in line:
                        splitline = line.split()
                        if id in genomes_cluster_dict.keys():
                            if splitline[1] not in genomes_cluster_dict[id]:
                                genomes_cluster_dict[id].append(splitline[1])
                        if id not in genomes_cluster_dict.keys():
                            genomes_cluster_dict[id].append(splitline[1])
    return genomes_cluster_dict



validchoice = False
while validchoice == False:
    datasetprompt = input("Contents page for both or comp genome dataset?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        mastloc = "Thesis_project_data\\Data\\Processed_data\\mastout"
        outputloc = "Thesis_project_data\\Data\\Processed_data\\IGR_association_analysis_both_genomes\\Masttxt_cluster_occurrence_pergenome.txt"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        mastloc =  "Thesis_project_data\\Data\\Processed_data\\mastout"
        outputloc = "Thesis_project_data\\Data\\Processed_data\\IGR_association_analysis_both_genomes\\Masttxt_cluster_occurrence_pergenome.txt"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")
genomeids = os.listdir(mastloc)
genomeclustdict = mastreader(mastloc,genomeids)
with open(outputloc, "w") as file:
    file.write("FIND WHICH MAST.TXT FILES CONTAIN WHICH IGR USING THIS\n")
    file.write("#######################################################\n")
    for genome in genomeclustdict.keys():
        print(f"Writing lines in output for {genome}")
        clusterline = '\t'.join(genomeclustdict[genome])
        file.write(f"{genome}\t{clusterline}\n")
