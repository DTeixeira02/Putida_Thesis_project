#Cluster-gene association comparator_V3
#CHANGES:
#For each TF it encodes, it writes a list of identified binding sites for that TF. For the TFs it doesnt APPEAR to encode, it still lists the binding sites but states TF is not present or unknown gene
#Given that multiple outputs occur for each genome, each genomes outputs are grouped into files listed under their name
#EXTRA FILE: Sectioned by TF, lists each identified clusterID that was binding site, and lists the percentages of the genes associated with those genes 
import os
import re
from collections import defaultdict
import pandas as pd
import time

############################################################
#                   DEFINING FUNCTIONS                     #
############################################################

def binding_site_lister(IGRcsv,hmmerloc): #READS HMMER BINDING SITE PREDICTIONS, TRANSLATES TO CLUSTERID
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
        print(f"Task at {(rowcounter/denominator)*100} % completion", end="\r")
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


def genelistreader(genelistdir):
    genome_to_namedgene = defaultdict(list)
    genome_to_IGR_op_genes = defaultdict(list) #Key = genomeid, values = operons in genome + IGR regulating it
    gene_to_IGRtagGENOME = defaultdict(list)
    genelist_list = os.listdir(genelistdir) #I am a master of naming variables
    print(f"Beginning Genelist readings ({len(genelist_list)} total)")
    counter = 0
    for genelist in genelist_list:
        print(f"Task at {100*counter/len(genelist_list)} % completion",end="\r")
        counter = counter+1
        beginreading = False
        id = re.findall("GENELIST_(.*)\.txt",genelist)[0]
        with open(str(genelistdir+"\\"+genelist), "r") as file:
            firstoperon = True
            for line in file:
                if beginreading == True:
                    splitline = line.split()
                    if len(splitline) == 4:
                        gene_to_IGRtagGENOME[splitline[2]].append([splitline[3],splitline[1],id]) #key = gene name, value = [clusterID,locustag,genomeID]
                        if firstoperon == False:
                            if currentoperon == splitline[0]:
                                opgenes.append(splitline[2])
                            if currentoperon != splitline[0]:
                                genome_to_IGR_op_genes[id].append([currentigr,opgenes])
                                currentoperon = splitline[0]
                                currentigr = splitline[3]
                                opgenes = [splitline[2]]
                        if firstoperon == True:
                            currentoperon = splitline[0]
                            currentigr = splitline[3]
                            opgenes = [splitline[2]]
                            firstoperon = False
                        if id in genome_to_namedgene.keys():
                            if splitline[2] not in genome_to_namedgene[id] and splitline[2] != "Unknown_Gene":
                                genome_to_namedgene[id].append(splitline[2])
                        if id not in genome_to_namedgene.keys():
                            if splitline[2] != "Unknown_Gene":
                                genome_to_namedgene[id].append(splitline[2])
                    else:
                        print(f"Potential error in {genelist} at line {line}")
                if "LOCUS TAG TO GENE LIST" in line:
                    beginreading = True
            genome_to_IGR_op_genes[id].append([currentigr,opgenes])
    print("Genelists read, moving to association comparisons")
    return gene_to_IGRtagGENOME, genome_to_namedgene, genome_to_IGR_op_genes

def IGR_to_gene_association_comparator(genome_to_IGRop):
    print("Beginning IGR to Operon association analysis")
    IGR_to_gene_assoc_data = defaultdict(list)
    IGR_to_dwnstrm_ops = defaultdict(list) 
    for genome in genome_to_IGRop.keys():
        for igr,operon in genome_to_IGRop[genome]:
            IGR_to_dwnstrm_ops[igr].append(operon)
    for igr in IGR_to_dwnstrm_ops.keys():
        distinct_op_occurrences = {} #Each time a gene occurs in a new operon, the value associated with that gene name gets incremented in this dict. Key = gene, Val = op occurrences
        single_gene_occ = [] #Records the name of every gene associated with the IGR only once, even if its assoc multiple times to same igr
        total_gene_occ = [] #Records every single occurrence of each gene assoc to current IGR
        total_assoc_to_op = len(IGR_to_dwnstrm_ops[igr])
        oplengths = []
        for operon in IGR_to_dwnstrm_ops[igr]:
            oplengths.append(len(operon))
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
                if gene in genes_checked_this_operon and gene != "Unknown_Gene":
                    print(f"{gene} occurs twice in the same operon associated to {igr}")
                genes_checked_this_operon.append(gene)
        avglength = (sum(oplengths)/len(oplengths))
        for occ in single_gene_occ:
            percent_assoc = 100 * distinct_op_occurrences[occ]/total_assoc_to_op
            IGR_to_gene_assoc_data[igr].append([occ,distinct_op_occurrences[occ],percent_assoc,total_assoc_to_op,avglength]) #total assoc to op only included so i dont have to pass the IGR_to_dwnstream_genes dict for len()
    print("Completed IGR to Operon aassociation analysis")
    return IGR_to_gene_assoc_data

def gene_to_IGR_association_comparator(genetoIGG): #IGG shorthand for IGR gene genome, variable is the dict produced by the above function
    print("Beginning Gene to IGR association analysis")
    significant_gene_freq_dict = defaultdict(list)
    gene_freq_dict = defaultdict(list)
    for gene in genetoIGG.keys():
        if gene != "Unknown_Gene":
            singleclusteroccs = [] #for each gene, holds all clusterIDs associated with it WITHOUT DUPLICATES - Used to tell .count() function what to count, without entering same clusterID twice
            totalclusteroccs = [] #for each gene, holds all clusterIDs associated with it WITH DUPLICATES - For clusterID counting purposes
            cluster_count = [] #nested list, each entry in format [cluster, numberoccurrencesofcluster,percentoftotaloccs]
            for cluster,tag,genome in genetoIGG[gene]:
                totalclusteroccs.append(cluster)
                if cluster not in singleclusteroccs:
                    singleclusteroccs.append(cluster)
            denominator = len(totalclusteroccs)
            expected_if_random_dist = ((denominator/len(singleclusteroccs))/denominator)*100
            for cluster in singleclusteroccs:
                numberassocs = totalclusteroccs.count(cluster)
                fractionassoc = (numberassocs/denominator)*100
                cluster_count.append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist]) 
                if fractionassoc > (15+expected_if_random_dist) and len(singleclusteroccs) != 1: #15 is arbitrary, checks if fraction assoc is significantly higher than expected.
                    significant_gene_freq_dict[gene].append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist]) 
                if len(singleclusteroccs) == 1 and denominator >= 50:
                    significant_gene_freq_dict[gene].append([cluster,denominator,numberassocs,fractionassoc,expected_if_random_dist])
            gene_freq_dict[gene].extend(cluster_count)
    print("Completed Gene to IGR association analysis")
    return significant_gene_freq_dict, gene_freq_dict

def binding_site_lister_pergenome(bs_to_tf, genome_to_IGRop, currentgenome, genome_to_namedgenes,conf_tf,conf_not_tf): #Runs once per genome - checks each tf binding IGR and divides them into ones the genome produces TFs for, and ones it doesnt. Also lists TFs the genome produces without apparent binding sites
    tf_list = ["AauR","AlgR","AmpR","AmrZ","Anr","ArgR","AtzR","BrlR","CatR","CifR","ColR","Cra","CueR","DesT","DNR","ExsA","FleQ","Fur","HexR","HrpL","IHF","IscR","LasR","LexA","MexT","MvfR","NarL","NtrC","OhrR","OxyR","PA2206","PcaR","PchR","PhhR","PhoB","PhzR","PltR","PsrA","PtxR","PtxS","PvdS","PyeR","RcsB","RhpR","RpoN","RsaL","TodT","TrpI","Vfr","VqsM","VqsR","Zur"]
    genometfs = []
    bsIGRs_with_TFs = {} #bs stands for binding site, dict keyed by IGR, val by TF that binds it. Contains every IGR in the genome that binds a TF that the genome encodes
    bsIGRs_without_TFs = defaultdict(list) #bs stands for binding site, dict keyed by IGR and TF that binds, contains every IGR that is a TF binding site BUT the genome DOESNT ENCODE THE TF
    TFs_without_IGRs = [] #Contains every TF the genome encodes from collecTF that doesnt appear to have a binding site in the genome
    remaining_IGRs = [] #Contains all IGRs that bind TFs we havent got data for(i.e all other IGRs not covered in the other dicts)
    genelist = []
    for namedgene in genome_to_namedgenes[currentgenome]:
        if namedgene not in genelist:
            genelist.append(namedgene)
        for tf in tf_list:
            if tf not in genometfs:
                if tf.lower() == namedgene.lower(): #Checks if the transcription factor name is identical to the current named gene, .lower() used to account for case sensitivity
                    genometfs.append(tf)
                    break #Named gene has been associated with a transcription factor, so continuing to try associate named gene with other TFs is unnecessary
                if tf.lower() in namedgene.lower() and namedgene.lower() in conf_tf: #If the current named gene is similar to the current TF and the current named gene is manually confirmed to be a TF
                    genometfs.append(tf)
                    break
                if tf.lower() in namedgene.lower() and namedgene.lower() in conf_not_tf:
                    break
                if tf.lower() in namedgene.lower() and tf.lower() != namedgene.lower() and namedgene.lower() not in conf_tf and namedgene.lower() not in conf_not_tf: #Essentially, if the named gene seems similar to a transcription factor and has yet to be confirmed as one or not, this if statement is executed
                    validchoice = False
                    while validchoice == False:
                        print("")
                        userprompt = input(f"{namedgene} has similar name to transcription factor {tf}, is {namedgene} a transcription factor?(Y/N):")
                        if userprompt == "Y" or userprompt == "y":
                            genometfs.append(tf)
                            conf_tf.append(namedgene.lower())
                            print(f"{namedgene} confirmed as isoform or subunit of {tf}")
                            validchoice = True
                        elif userprompt == "N" or userprompt == "n":
                            conf_not_tf.append(namedgene.lower())
                            validchoice = True
                    break
    for igr,operon in genome_to_IGRop[currentgenome]:
        if igr in bs_to_tf.keys():
            if igr not in bsIGRs_with_TFs.keys() and igr not in bsIGRs_without_TFs.keys():
                tffound = False
                for tf in bs_to_tf[igr]:
                    if tf in genometfs:
                        if tffound == True:
                            print(f"For {igr} in {currentgenome}, multiple TFs can bind and BOTH are encoded in this genome {bs_to_tf[igr]}")
                        tffound = True
                if tffound == True:
                    bsIGRs_with_TFs[igr] = bs_to_tf[igr]
                elif tffound == False:
                    bsIGRs_without_TFs[igr] = bs_to_tf[igr]
            elif igr in bsIGRs_with_TFs:
                print(f"{igr} in {currentgenome} has already been associated to {bsIGRs_with_TFs[igr]}; Attempted association to {tf}")
        if igr not in bs_to_tf.keys():
            remaining_IGRs.append(igr)
    for tf in genometfs:
        if tf not in bsIGRs_with_TFs.values():
            TFs_without_IGRs.append(tf)
    return bsIGRs_with_TFs,bsIGRs_without_TFs,TFs_without_IGRs,remaining_IGRs,conf_tf,conf_not_tf

############################################################
#                   USER INPUT                             #
############################################################
validchoice = False
while validchoice == False:
    datasetprompt = input("What dataset are you wanting to analyse IGR/Gene associations for?('both' or 'comp'):")
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
#                 MAIN FUNCTION EXECUTION                  #
############################################################
bs_to_tf_dict = binding_site_lister(piggy_igr_loc,hmmerloc)
gene_to_igrgenegenome_dict,genome_to_gene,genome_to_op_dict = genelistreader(genelistloc)
sig_gene_dict, gene_dict = gene_to_IGR_association_comparator(gene_to_igrgenegenome_dict)
igr_assoc_data_dict = IGR_to_gene_association_comparator(genome_to_op_dict)

############################################################
#                  FILE WRITE EXECUTION                    #
############################################################

unformattedids = os.listdir(genelistloc)
formattedids = []
progressdenom = len(unformattedids)
progressnumer = 0
for ufid in unformattedids:
    id = re.findall("GENELIST_(.*)\.txt",ufid)[0]
    formattedids.append(id)
    foldertomake = str(outfilestem+"\\"+id)
    try:
        os.mkdir(foldertomake)
    except FileExistsError:
        print(f"File for {id} already exists in {outfilestem}, contents will be overwritten")
conf_tf = []
conf_not_tf = []
for id in formattedids:
    bsIGRs_with_TFs,bsIGRs_without_TFs,TFs_without_IGRs,remaining_IGRs,conf_tf,conf_not_tf = binding_site_lister_pergenome(bs_to_tf_dict, genome_to_op_dict, id, genome_to_gene,conf_tf,conf_not_tf)
    progressnumer = progressnumer + 1
    outputloc = str(outfilestem+"\\"+id)
    print(f"Writing association analysis for {id}; Overall progress {round(((progressnumer/progressdenom)*100),2)}", end="\r")
    with open(str(outputloc+"\\Gene_to_IGR_assoc_analysis.txt"),"w") as file:
        if id in genome_to_gene.keys():
            file.write(f"{id} GENE TO IGR ASSOCIATIONS\n\n#######################################################################\n")
            file.write("SIGNIFICANT GENE CLUSTER ASSOCIATIONS OCCURRING >100x WITH FRACTION ASSOCIATED +15% HIGHER THAN EXPECTED\n")
            file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
            for gene in genome_to_gene[id]:
                clusteroccs = []
                for igr,tag,genome in gene_to_igrgenegenome_dict[gene]:
                    if genome == id:
                        clusteroccs.append(igr)
                for cluster,denominator,numberassoc,fracassoc,expected_dist in sig_gene_dict[gene]:
                    if cluster in clusteroccs:
                        file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{round(fracassoc,3)}\t{round(expected_dist,3)}\n")
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
                        file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{round(fracassoc,3)}\t{round(expected_dist,3)}\n")
        elif id not in genome_to_gene.keys():
            print(f"An error has occurred, {id} doesn't have any listed named genes")
    with open(str(outputloc+"\\IGR_TF_binding_analysis.txt"),"w") as file:
        file.write(f"{id} IGR TF BINDING ANALYSIS\n\n#######################################################################\n")
        file.write(f"TF-Binding IGRs where genome encodes TF\n")
        file.write("IGR\tTF-BOUND\n")
        for igr in bsIGRs_with_TFs.keys():
            if len(bsIGRs_with_TFs[igr]) > 1:
                tf = ",".join(bsIGRs_with_TFs[igr])
            elif len(bsIGRs_with_TFs[igr]) == 1:
                tf = bsIGRs_with_TFs[igr][0]
            file.write(f"{igr}\t{tf}\n")
        file.write("\n#######################################################################\n")
        file.write(f"TF-Binding IGRs where genome doesn't encode TF\n")
        file.write("IGR\tTF-BOUND\n")
        for igr in bsIGRs_without_TFs:
            if len(bsIGRs_without_TFs[igr]) > 1:
                tf = ",".join(bsIGRs_without_TFs[igr])
            elif len(bsIGRs_without_TFs[igr]) == 1:
                tf = bsIGRs_without_TFs[igr][0]
            file.write(f"{igr}\t{tf}\n")
        file.write("\n#######################################################################\n")
        file.write(f"TFs encoded by genome lacking binding site\n")
        file.write("TRANSCRIPTION FACTOR\n")
        for tf in TFs_without_IGRs:
            file.write(f"{tf}\n")
        file.write("\n#######################################################################\n")
        file.write(f"Remaining IGRs with unknown TF binding profiles\n")
        file.write("IGR\tOCC_COUNT\n")
        remainingappended = []
        for igr in remaining_IGRs:
            if igr != "NO_CLUSTER_ASSOCIATION_FOUND" and igr not in remainingappended:
                igrcount = remaining_IGRs.count(igr)
                file.write(f"{igr}\t{igrcount}\n")
                remainingappended.append(igr)
        
print("\nWriting the collated file summarising all gene to IGR associations across all genomes")
with open(str(outfilestem+"\\Collated_gene_to_igr_associations.txt"), "w") as file:
    file.write("SIGNIFICANT GENE CLUSTER ASSOCIATIONS OCCURRING >100x WITH FRACTION ASSOCIATED +15% HIGHER THAN EXPECTED\n")
    file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
    for gene in sig_gene_dict.keys():
        for cluster,denominator,numberassoc,fracassoc,expected_dist in sig_gene_dict[gene]:
            if numberassoc > 100:
                file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{round(fracassoc,3)}\t{round(expected_dist,3)}\n")
    file.write("\n#######################################################################\n")
    file.write("ALL GENE TO IGR ASSOCIATIONS ACROSS ALL GENOMES\n")
    file.write("GENE\tIGR\tGENE_COUNT\tASSOC_COUNT\tASSOC_TOTAL(%)\tASSOC_EXPECTED(%)\n")
    for gene in gene_dict:
        for cluster,denominator,numberassoc,fracassoc,expected_dist in gene_dict[gene]:
            file.write(f"{gene}\t{cluster}\t{denominator}\t{numberassoc}\t{round(fracassoc,3)}\t{round(expected_dist,3)}\n")

print("Writing collated IGR to gene associations from across all genomes")
with open(str(outfilestem+"\\Collated_bsigr_to_gene_associations.txt"),"w") as file:
    file.write("This file lists every transcription factor binding IGR proposed to regulate an operon, and statistics for the operon composition for each regulatory occurrence of the IGR\n")
    file.write("If IGR occurrence isn't associated to an operon, it isnt included when calculating statistics\n")
    file.write("IGR = ClusterID, AVG_OP_LEN = Average number of genes per operon to 3 D.P,TF = Name of TF(s) IGR binds separated by comma, TOTAL_ASSOC = Number of IGR occurrences associated to operon(s), GENE = Name of gene if known, OCC_IN_OPERON = Number of operons associated to IGR bearing gene\n")
    file.write("IGR\tAVG_OP_LEN\tTF\tTOTAL_ASSOCS\tGENE\tOCC_IN_OPERON\tPERCENT_OPERONS\n")
    for igr in igr_assoc_data_dict.keys():
        if igr in bs_to_tf_dict.keys():
            if len(bs_to_tf_dict[igr]) > 1:
                currentTF = ','.join(bs_to_tf_dict[igr])
            elif len(bs_to_tf_dict[igr]) == 1:
                currentTF = bs_to_tf_dict[igr][0]
            for gene,occ_count,percent_ops,total_assocs,avglen in igr_assoc_data_dict[igr]:
                file.write(f"{igr}\t{round(avglen,3)}\t{currentTF}\t{total_assocs}\t{gene}\t{occ_count}\t{round(percent_ops,3)}\n")

