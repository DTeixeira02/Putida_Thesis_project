#Operon to IGR edge creator VERSION3
#Reads .gbk files for each genome to ID gene regulators that we have found binding sites for from collecTF
#IDs which clusters are likely binding sites for a given TF based on binding site sequences from collecTF processed by hmmer
#Forms edges between operons expressing a given TF and the IGR(s) in the same genome that TF likely binds to
#Note - doesn't account for all TFs in the genome, only those we have info on from collecTF
#Note - IGRs predicted as TF binding sites aren't necessarily accurate

#VERSION 2 NOTES
#Previously, isoforms of transcription factors from the .gbk were ignored such as cueR_1/cueR_2/cueR_3... only cueR was recognised
#This version attempts to account for as many isoforms as possible

#VERSION 3 NOTES
#Applies same changes as V2
#If IGR found as TF binding site in one genome, all occurrences of that IGR are also assumed binding sites across all genomes

import networkx as nx
from collections import defaultdict
import os
import re
import pandas as pd



def genelistreader(currentgenelist, conf_tf_list, conf_not_tf_list):
    tf_list = ["AauR","AlgR","AmpR","AmrZ","Anr","ArgR","AtzR","BrlR","CatR","CifR","ColR","Cra","CueR","DesT","DNR","ExsA","FleQ","Fur","HexR","HrpL","IHF","IscR","LasR","LexA","MexT","MvfR","NarL","NtrC","OhrR","OxyR","PA2206","PcaR","PchR","PhhR","PhoB","PhzR","PltR","PsrA","PtxR","PtxS","PvdS","PyeR","RcsB","RhpR","RpoN","RsaL","TodT","TrpI","Vfr","VqsM","VqsR","Zur"]
    tag_to_op_dict = {}
    tag_to_igr_dict = {}
    tf_to_tag_dict = defaultdict(list)
    namedgene_tags_dict = defaultdict(list)
    beginreading = False
    print(f"Beginning read for {currentgenelist}")
    with open(currentgenelist, "r") as file:
        for line in file:
            if beginreading == True:
                splitline = line.split()
                if len(splitline) == 4:
                    tag_to_op_dict[splitline[1]] = splitline[0]
                    tag_to_igr_dict[splitline[1]] = splitline[3]
                    if splitline[2] != "Unknown_Gene":
                        namedgene_tags_dict[splitline[2]].append(splitline[1])
                else:
                    print(f"Potential error in {currentgenelist} at line {line}")
            if "LOCUS TAG TO GENE LIST" in line:
                beginreading = True
    for named_gene in namedgene_tags_dict.keys(): #dict keyed by gene names, values are lists of locus tags bearing that gene name {genename:[loc1,loc2...]}
        for tf in tf_list: #Loops for every Transcription factor I got binding data for from collecTF
            if tf.lower() == named_gene.lower(): #Checks if the transcription factor name is identical to the current named gene, .lower() used to account for case sensitivity
                tf_to_tag_dict[tf].extend(namedgene_tags_dict[named_gene])
                break #Named gene has been associated with a transcription factor, so continuing to try associate named gene with other TFs is unnecessary
            if tf.lower() in named_gene.lower() and named_gene.lower() in conf_tf_list: #If the current named gene is similar to the current TF and the current named gene is manually confirmed to be a TF
                tf_to_tag_dict[tf].extend(namedgene_tags_dict[named_gene])
                break
            if tf.lower() in named_gene.lower() and named_gene.lower() in confirmed_not_tf_list:
                break
            if tf.lower() in named_gene.lower() and tf.lower() != named_gene.lower() and named_gene.lower() not in conf_tf_list and named_gene.lower() not in conf_not_tf_list: #Essentially, if the named gene seems similar to a transcription factor and has yet to be confirmed as one or not, this if statement is executed
                validchoice = False
                while validchoice == False:
                    userprompt = input(f"{named_gene} has similar name to transcription factor {tf}, is {named_gene} a transcription factor?(Y/N):")
                    if userprompt == "Y" or userprompt == "y":
                        tf_to_tag_dict[tf].extend(namedgene_tags_dict[named_gene])
                        conf_tf_list.append(named_gene.lower())
                        print(f"Appended {named_gene} data under {tf} in tf_to_tag_dict")
                        validchoice = True
                    elif userprompt == "N" or userprompt == "n":
                        confirmed_not_tf_list.append(named_gene.lower())
                        validchoice = True
                break
    print(f"{currentgenelist} read, moving onto graph updating")
    return tag_to_op_dict, tag_to_igr_dict, tf_to_tag_dict

def binding_site_lister(IGRcsv,hmmerloc): #READS HMMER BINDING SITE PREDICTIONS, TRANSLATES TO CLUSTERID
    tf_to_bs_dict = defaultdict(list) #Dictionary where keys are TFs and values are the IGRs they bind
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
                    if currentTF in tf_to_bs_dict.keys():
                        if Locus_string_to_IGR[splitline[0]] not in tf_to_bs_dict[currentTF]:
                            tf_to_bs_dict[currentTF].append(Locus_string_to_IGR[splitline[0]])
                    elif currentTF not in tf_to_bs_dict.keys():
                        tf_to_bs_dict[currentTF].append(Locus_string_to_IGR[splitline[0]])       
    return tf_to_bs_dict


##########################
##########################
##########BODY############
##########################
##########################
validchoice = False
while validchoice == False:
    datasetprompt = input("What dataset are you working with?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        genelistloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_both_genomes"
        graphfolder = "Internship_data\\Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_bothgenomes"
        outputgraphfolder = "Thesis_project_data\\Data\\Processed_data\\V3_Graphml_full_regulatory_network_bothgenomes\\Regulatory_network_graph_"
        igrcsvloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        genelistloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_comp_genomes"
        graphfolder = "Internship_data\\Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_compgenomes"
        outputgraphfolder = "Thesis_project_data\\Data\\Processed_data\\V3_Graphml_full_regulatory_network_compgenomes\\Regulatory_network_graph_"
        igrcsvloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")
igrgeneloc = "Internship_data\\Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt"
hmmertabloc = "Internship_data\\Data\\Hmmer_data\\nhmmerout"
unformattedids = os.listdir(genelistloc)

confirmed_tf_list = []
confirmed_not_tf_list = []
TF_to_bindsite_dict = binding_site_lister(igrcsvloc,hmmertabloc)
for x in unformattedids:
    id = re.findall("GENELIST_(.*).txt",x)[0]
    print(f"Beginning graph completion process for {id}")
    tfs_to_encoding_operons_dict = defaultdict(list)
    tag_op_dict, tag_igr_dict, tf_tags_dict = genelistreader(str(genelistloc+"\\"+x),confirmed_tf_list,confirmed_not_tf_list)
    for tf in tf_tags_dict.keys():
        for tag in tf_tags_dict[tf]:
            if tag_op_dict[tag] not in tfs_to_encoding_operons_dict[tf]:
                try:
                    tfs_to_encoding_operons_dict[tf].append(tag_op_dict[tag])
                except KeyError: #Incase the gbkreader somehow returned a locustag that doesnt exist in genome 'x'
                    if tag not in tag_op_dict.keys():
                        print(f"ERROR - {tag} is not a locus tag found in {x}")
    originalgraph = nx.read_graphml(str(graphfolder+"\\IGR_OPERON_GML-"+id+".graphml"))
    currentgraphnodes = list(originalgraph.nodes())
    for tf in tfs_to_encoding_operons_dict.keys(): #For transcription factor from collecTF data that are known to be encoded by operon in current genome 
        if tf in TF_to_bindsite_dict.keys(): #If that TF has a binding site in any genome(assumes IGR that binds tf in one genome binds same TF in all genomes)
            for operon in tfs_to_encoding_operons_dict[tf]: #For every operon encoding that TF
                for igr in TF_to_bindsite_dict[tf]: #For every IGR identified as a binding site in all genomes
                    if operon in currentgraphnodes and igr in currentgraphnodes: #Makes sure both nodes exist, given that igr data is coming from all genomes and not just currentgenome
                        originalgraph.add_edge(operon, igr) #Forms edge from aforementioned operon(s) to aforementioned IGR(s)
                    elif operon not in currentgraphnodes:
                        print(f"ERROR - {operon} doesn't exist as node in graph for {id}")
        elif tf not in TF_to_bindsite_dict.keys(): #If TF doesnt have a binding site in any genome
            print(f"{tf} binding site not found in any genome")
    for node in list(originalgraph.nodes()):
        if "Cluster_" in node:
            originalgraph.nodes[node]["color"] = "#33D1FF" #sets node colour for IGRs to a light blue
        elif "Operon" in node:
            originalgraph.nodes[node]["color"] = "#FFE033" 
    nx.write_graphml(originalgraph, str(outputgraphfolder+id+".graphml"))
