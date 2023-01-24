#IGR-Gene_association_variation_analysis
import os
from collections import defaultdict
############################################################
#                   DEFINING FUNCTIONS                     #
############################################################
def gene_igr_assoc_reader(col_gene_igr_txt):
    gene_to_igrstats = defaultdict(list)
    with open(col_gene_igr_txt, "r") as file:
        begin_reading = False
        for line in file:
            if "ALL GENE TO IGR ASSOCIATIONS ACROSS ALL GENOMES" in line:
                begin_reading = True
            if begin_reading == True:
                splitline = line.split()
                if "Cluster_" in splitline[1]:
                    gene_to_igrstats[splitline[0]].append([splitline[1],splitline[2],splitline[3],splitline[4],splitline[5]]) #GENE:[IGR,GENE_COUNT,ASSOC_COUNT,ASSOC_TOTAL(%),ASSOC_EXPECTED(%)]
    return gene_to_igrstats            


def igr_gene_assoc_reader(col_igr_gene_txt,gene_to_igrstats):
    bsigr_to_tf_dict = {}
    bsigr_to_opgenes = defaultdict(list)
    gene_to_bsigrs = defaultdict(list)
    genes_with_multiple_bsigrs = defaultdict(list)
    with open(col_igr_gene_txt, "r") as file:
        first_cluster = True
        for line in file:
            splitline = line.split()
            if len(splitline) == 7 and "Cluster_" in splitline[0]:
                if first_cluster == True:
                    currentcluster = splitline[0]
                    first_cluster = False
                if splitline[0] == currentcluster:
                    bsigr_to_opgenes[currentcluster].append(splitline[4])
                if splitline[0] != currentcluster:
                    currentcluster = splitline[0]
                    bsigr_to_opgenes[currentcluster].append(splitline[4])
                if splitline[0] not in bsigr_to_tf_dict.keys():
                    bsigr_to_tf_dict[splitline[0]] = splitline[2]
                if splitline[0] in bsigr_to_tf_dict.keys():
                    if splitline[2] not in bsigr_to_tf_dict[splitline[0]]:
                        print(f"Potential error with {splitline[0]}: Different instances bind different TFs, {bsigr_to_tf_dict[splitline[0]]} vs {splitline[2]}")
                if splitline[4] in gene_to_bsigrs.keys():
                    if splitline[0] not in gene_to_bsigrs[splitline[4]]:
                        gene_to_bsigrs[splitline[4]].append(splitline[0])
                if splitline[4] not in gene_to_bsigrs.keys():
                    gene_to_bsigrs[splitline[4]].append(splitline[0])
    for gene in gene_to_bsigrs: #For every gene found associated to an IGR for which TF binding profile known
        if len(gene_to_bsigrs[gene]) > 1: #If the gene is found associated to multiple different IGRs
            tmptfs = []
            tmpdata = []
            for igr in gene_to_bsigrs[gene]: #For each IGR the gene is found assoc to
                if bsigr_to_tf_dict[igr] not in tmptfs: #If the TF that binds the current IGR is different to any past TFs binding previous IGRs associated to this gene
                    for clusterid,genecount,assoccount,percentassoc,assocexp in gene_to_igrstats[gene]: #Getting stats across genomes for this gene to IGR association
                        if clusterid == igr:
                            tmptfs.append(bsigr_to_tf_dict[igr])
                            tmpdata.append([igr,bsigr_to_tf_dict[igr],genecount,percentassoc,assocexp]) #[ClustID, TF-bound, Number_times_gene_occs, Percent_time_gene_occs_with_ClustID, Percent_times_gene_expected_assoc_ClustID]
            if len(tmptfs)>1: #If after all of the above, the gene was found to be regulated by more than one TF across all genomes
                genes_with_multiple_bsigrs[gene].extend(tmpdata) #Gene key, val = Each cluster it assoc to and the stats for gene to IGR assoc
    return genes_with_multiple_bsigrs

############################################################
#                   USER INPUT                             #
############################################################
validchoice = False
while validchoice == False:
    datasetprompt = input("What dataset are you wanting to analyse association variation for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        outfilestem = "Thesis_project_data\\Data\\Processed_data\\IGR_assoc_analysis\\IGR_assoc_analysis_bothgenomes"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        outfilestem = "Thesis_project_data\\Data\\Processed_data\\IGR_assoc_analysis\\IGR_assoc_analysis_compgenomes"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")
############################################################
#                 MAIN FUNCTION EXECUTION                  #
############################################################
print("Beginning collated Gene to IGR association file read")
gene_igr_stats_dict = gene_igr_assoc_reader(str(outfilestem+"\\Collated_gene_to_igr_associations.txt"))
print("Beginning collated Transcription-factor binding IGR to Gene association file read")
gene_multi_bsigrs_dict = igr_gene_assoc_reader(str(outfilestem+"\\Collated_bsigr_to_gene_associations.txt"),gene_igr_stats_dict)

############################################################
#                  FILE WRITE EXECUTION                    #
############################################################
print("Writing results to file")
with open(str(outfilestem+"\\Gene_to_bsigr_variation_analysis.txt"),"w") as file:
    file.write("LISTS ALL GENES FROM Collated_bsigr_to_gene_associations.txt THAT ARE FOUND ASSOCIATED TO 2 OR MORE IGRS THAT BIND DIFFERENT TRANSCRIPTION FACTORS\nFORMATTED ACROSS LINES: GENE NAME AS HEADER, FOLLOWING LINES AND STATISTICS FOR IGRS GENE FOUND ASSOCIATED TO\nIGR STATS FORMAT: CLUSTERID, TF-BOUND, GENE_COUNT_ACROSS_GENOMES, PERCENTAGE_GENE_COUNT_ASSOC_TO_IGR, EXPECTED_PERCENATGE\n")
    file.write("###############################################################################################################################\n")
    for gene in gene_multi_bsigrs_dict.keys():
        file.write(f"\n{gene}\n")
        for igr,tf,numoccs,percentassoc,assocexp in gene_multi_bsigrs_dict[gene]:
            file.write(f"{igr}\t{tf}\t{numoccs}\t{percentassoc}\t{assocexp}\n")
print("Task complete")