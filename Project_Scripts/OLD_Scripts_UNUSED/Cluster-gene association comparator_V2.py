#Cluster-gene association comparator_V2
#CHANGES:
#For each TF it encodes, it writes a list of identified binding sites for that TF. For the TFs it doesnt APPEAR to encode, it still lists the binding sites but states TF is not present or unknown gene
#Given that multiple outputs occur for each genome, each genomes outputs are grouped into files listed under their name
#EXTRA FILE: Sectioned by TF, lists each identified clusterID that was binding site, and lists the percentages of the genes associated with those genes 
import os
import re
from collections import defaultdict
import pandas as pd
import time

def binding_site_lister(IGRcsv,hmmerloc):
    bs_to_tf_dict = defaultdict(list) #Dictionary where keys are IGRs and values are the TFs they bind
    Locus_string_to_IGR = {}
    rowcounter = 0
    csvdata = pd.read_csv(IGRcsv, low_memory=False)
    clusterIDs = csvdata.iloc[:,0]
    denominator = len(clusterIDs)
    print("Beginning IGR csv reading process")
    for cluster in clusterIDs:
        rowdata = csvdata.iloc[rowcounter,:]
        for cell in rowdata:
            if pd.isnull(cell) == False:
                if "GCA_" in str(cell):
                    splittingdoubles = cell.split()
                    for item in splittingdoubles:
                        if item in Locus_string_to_IGR.keys():
                            print(f"Error occurred, {item} already associated with {Locus_string_to_IGR[item]}, cannot associate to {cluster}")
                        if item not in Locus_string_to_IGR.keys():
                            Locus_string_to_IGR[item] = cluster
        rowcounter = rowcounter+1
        print(f"Task at {(rowcounter/denominator)*100} % completion")
    filelist = os.listdir(hmmerloc)
    for hmmerpred in filelist:
        currentTF = re.findall("(.*)_pred.tab",hmmerpred)[0]
        with open(str(hmmerloc+"\\"+hmmerpred),"r") as file:
            for line in file:
                if "GCA_" in line:
                    splitline = line.split()
                    if Locus_string_to_IGR[splitline[0]] in bs_to_tf_dict.keys():
                        if currentTF not in bs_to_tf_dict[Locus_string_to_IGR[splitline[0]]]:
                            print(f"Potential error - {Locus_string_to_IGR[splitline[0]]} is a binding site for {bs_to_tf_dict[Locus_string_to_IGR[splitline[0]]]} AND APPARENTLY {currentTF}")
                            bs_to_tf_dict[Locus_string_to_IGR[splitline[0]]].append(currentTF)
                    elif Locus_string_to_IGR[splitline[0]] not in bs_to_tf_dict.keys():
                        bs_to_tf_dict[Locus_string_to_IGR[splitline[0]]].append(currentTF)             
    return bs_to_tf_dict


def genelistreader_igrtogene_assoc_comparator(genelistdir,bs_to_tf):
    genome_to_namedgene = defaultdict(list)
    gene_to_IGRtagGENOME = defaultdict(list)
    IGR_to_dwnstrm_genes = defaultdict(list)
    IGR_to_gene_assoc_data = defaultdict(list)
    genelist_list = os.listdir(genelistdir) #I am a master of naming variables
    print(f"Beginning Genelist readings ({len(genelist_list)} total)")
    for genelist in genelist_list:
        beginreading = False
        id = re.findall("GENELIST_(.*)\.txt",genelist)[0]
        with open(str(genelistdir+"\\"+genelist), "r") as file:
            tmpoperon = []
            firstoperon = True
            for line in file:
                if beginreading == True:
                    splitline = line.split()
                    if len(splitline) == 4:
                        gene_to_IGRtagGENOME[splitline[2]].append([splitline[3],splitline[1],id]) #key = gene name, value = [clusterID,locustag,genomeID]
                        if id in genome_to_namedgene.keys():
                            if splitline[2] not in genome_to_namedgene[id] and splitline[2] != "Unknown_Gene":
                                genome_to_namedgene[id].append(splitline[2])
                        if id not in genome_to_namedgene.keys():
                            if splitline[2] != "Unknown_Gene":
                                genome_to_namedgene[id].append(splitline[2])
                        if splitline[3] in bs_to_tf.keys():
                            if firstoperon == False:
                                if splitline[0] == currentoperon:
                                    tmpoperon.append(splitline[2])
                                elif splitline[0] != currentoperon:
                                    IGR_to_dwnstrm_genes[currentcluster].append(tmpoperon)
                                    tmpoperon = [splitline[2]]
                                    currentcluster = splitline[3]
                                    currentoperon = splitline[0]
                            elif firstoperon == True:
                                currentcluster = splitline[3]
                                tmpoperon.append(splitline[2])
                                currentoperon = splitline[0]
                                firstoperon = False      
                    else:
                        print(f"Potential error in {genelist} at line {line}")
                if "LOCUS TAG TO GENE LIST" in line:
                    beginreading = True
            if len(tmpoperon) != 0: #Previous code only appends operon once next operon assoc with binding site is found - Last of such operons in each file would be missed if not appended here
                IGR_to_dwnstrm_genes[currentcluster].append(tmpoperon) 
    print("Genelists read, moving to association comparisons")

    for igr in IGR_to_dwnstrm_genes.keys():
        distinct_op_occurrences = {} #Each time a gene occurs in a new operon, the value associated with that gene name gets incremented in this dict. Key = gene, Val = op occurrences
        single_gene_occ = [] #Records the name of every gene associated with the IGR only once, even if its assoc multiple times to same igr
        total_gene_occ = [] #Records every single occurrence of each gene assoc to current IGR
        total_assoc_to_op = len(IGR_to_dwnstrm_genes[igr])
        for operon in IGR_to_dwnstrm_genes[igr]:
            genes_checked_this_operon = [] #tmp list to ensure if same gene occurs twice in same operon, distinct_op_occurrences isnt incremented twice
            for gene in operon:
                total_gene_occ.append(gene)
                if gene not in single_gene_occ:
                    single_gene_occ.append(gene)
                if gene not in genes_checked_this_operon:
                    if gene in distinct_op_occurrences.keys():
                        distinct_op_occurrences[gene] = distinct_op_occurrences[gene] + 1
                    if gene not in distinct_op_occurrences.keys():
                        distinct_op_occurrences[gene] = 1
                if gene in genes_checked_this_operon:
                    print(f"{gene} occurs twice in the same operon associated to {igr}")
                genes_checked_this_operon.append(gene)
        for occ in single_gene_occ:
            percent_assoc = 100 * distinct_op_occurrences[occ]/total_assoc_to_op
            IGR_to_gene_assoc_data[igr].append([occ,distinct_op_occurrences[occ],percent_assoc,total_assoc_to_op]) #total assoc to op only included so i dont have to pass the IGR_to_dwnstream_genes dict for len()
    return gene_to_IGRtagGENOME, genome_to_namedgene, IGR_to_gene_assoc_data

def gene_to_IGR_association_comparator(genetoIGG): #IGG shorthand for IGR gene genome, variable is the dict produced by the above function
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
                    #print(f"{cluster} found associated with {numberassocs} copies of {gene} which occurred {denominator} times across all genomes")
                    significant_gene_freq_dict[gene].append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist]) 
                if len(singleclusteroccs) == 1 and denominator >= 50:
                    significant_gene_freq_dict[gene].append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist])
                    #print(f"{cluster} occurs {numberassocs} times and IS THE ONLY IGR ASSOCIATED WITH {gene} which occurred {denominator} times across all genomes")
            gene_freq_dict[gene].extend(cluster_count_list)
    return significant_gene_freq_dict, gene_freq_dict


############################################################
#                   USER INPUT                             #
############################################################
validchoice = False
while validchoice == False:
    datasetprompt = input("What dataset are you wanting to list genes for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        genelistloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_both_genomes"
        outfilestem = "Thesis_project_data\\Data\\Processed_data\\IGR_assoc_analysis\\IGR_assoc_analysis_bothgenomes"
        piggy_igr_loc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        hmmerloc = "Internship_data\\Data\\Hmmer_data\\nhmmerout"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        genelistloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_comp_genomes"
        outfilestem = "Thesis_project_data\\Data\\Processed_data\\IGR_assoc_analysis\\IGR_assoc_analysis_compgenomes"
        piggy_igr_loc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        hmmerloc = "Internship_data\\Data\\Hmmer_data\\nhmmerout"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")

############################################################
#                       BODY                               #
############################################################


unformattedids = os.listdir(genelistloc)
formattedids = []
bs_to_tf_dict = binding_site_lister(piggy_igr_loc,hmmerloc)
gene_to_igrgenegenome_dict,genome_to_gene, igr_assoc_data_dict = genelistreader_igrtogene_assoc_comparator(genelistloc,bs_to_tf_dict)
sig_gene_dict, gene_dict = gene_to_IGR_association_comparator(gene_to_igrgenegenome_dict)



for ufid in unformattedids:
    id = re.findall("GENELIST_(.*)\.txt",ufid)[0]
    formattedids.append(id)
    foldertomake = str(outfilestem+"\\"+id)
    try:
        os.mkdir(foldertomake)
    except FileExistsError:
        print(f"File for {id} already exists in {outfilestem}, contents will be overwritten")
for id in formattedids:
    outputloc = str(outfilestem+"\\"+id)
    time.sleep(3)
    print(f"Writing gene-igr associations for {id}")
    with open(str(outputloc+"\\GTI_assoc_analysis.txt"),"w") as file:
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

print("Writing the collated file summarising all gene to IGR associations across all genomes")
with open(str(outfilestem+"\\Collated_gene_to_igr_associations.txt"), "w") as file:
    file.write("SIGNIFICANT GENE CLUSTER ASSOCIATIONS OCCURRING >100x WITH FRACTION ASSOCIATED +15% HIGHER THAN EXPECTED\n")
    file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
    for gene in sig_gene_dict.keys():
        for cluster,denominator,numberassoc,fracassoc,expected_dist in sig_gene_dict[gene]:
            if numberassoc > 100:
                file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{fracassoc}\t{expected_dist}\n")
                #print(f"Look into {gene} downstream of {cluster} a total of {numberassoc} times ({fracassoc} of all occurrences, expected {expected_dist})")
    file.write("\n#######################################################################\n")
    file.write("ALL GENE TO IGR ASSOCIATIONS ACROSS ALL GENOMES\n")
    file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
    for gene in gene_dict:
        for cluster,denominator,numberassoc,fracassoc,expected_dist in gene_dict[gene]:
            file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{fracassoc}\t{expected_dist}\n")

print("Writing collated IGR to gene associations from across all genomes")
with open(str(outfilestem+"\\Collated_igr_to_gene_associations.txt"),"w") as file:
    file.write("This file lists every transcription factor binding IGR proposed to regulate an operon, an statistics for the operon composition for each regulatory occurrence of the IGR\n")
    file.write("If IGR occurrence isn't associated to an operon, it isnt included when calculating statistics\n")
    file.write("IGR = ClusterID, TF = Name of TF(s) IGR binds separated by comma, TOTAL_ASSOC = Number of IGR occurrences associated to operon(s), GENE = Name of gene if known, OCC_IN_OPERON = Number of operons associated to IGR bearing gene\n")
    file.write("IGR\tTF\tTOTAL_ASSOCS\tGENE\tOCC_IN_OPERON\tPERCENT_OPERONS\n")
    for igr in igr_assoc_data_dict.keys():
        if len(bs_to_tf_dict[igr]) > 1:
            currentTF = ','.join(bs_to_tf_dict[igr])
        elif len(bs_to_tf_dict[igr]) == 1:
            currentTF = bs_to_tf_dict[igr][0]
        for gene,occ_count,percent_ops,total_assocs in igr_assoc_data_dict[igr]:
            file.write(f"{igr}\t{currentTF}\t{total_assocs}\t{gene}\t{occ_count}\t{percent_ops}\n")