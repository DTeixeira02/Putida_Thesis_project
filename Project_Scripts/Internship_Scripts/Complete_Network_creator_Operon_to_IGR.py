#Operon to IGR edge creator
#Reads .gbk files for each genome to ID gene regulators that we have found binding sites for from collecTF
#IDs which clusters are likely binding sites for a given TF based on binding site sequences from collecTF processed by hmmer
#Forms edges between operons expressing a given TF and the IGR(s) in the same genome that TF likely binds to
#Note - doesn't account for all TFs in the genome, only those we have info on from collecTF
#Note - IGRs predicted as TF binding sites aren't necessarily accurate

import networkx as nx
from collections import defaultdict
import os
import re

def tabreader(tablist, tabloc): #loops over every hmmer .tab output, records location of igrs for igrs that tfs bind to in list of tuples [(genome, loc_tag, tf)], later find clusterID using loctag
    genelist = []
    for item in tablist:
        currentTF = re.findall("(.*)_pred.tab", item)[0]
        tab = str(tabloc + "\\" + item)
        with open(tab, "r") as file:
            for line in file.readlines():
                if "GCA_" in line:
                    currentline = line.split()
                    tags = currentline[0].split("_+_+_")
                    if tags[3] == "CO_R":
                        genetuple = (tags[0],tags[1],currentTF)
                        genelist.append(genetuple)
                    elif tags[3] == "CO_F":
                        genetuple = (tags[0],tags[2],currentTF)
                        genelist.append(genetuple)
                    elif tags[3] == "DP":
                        genetuple = (tags[0],tags[1],currentTF)
                        genelist.append(genetuple)
                        genetuple = (tags[0],tags[2],currentTF)
                        genelist.append(genetuple)
                    elif tags[3] == "DT":
                        genetuple = (tags[0],"N/A-IGRin/isdoubleterminator",currentTF)
                        genelist.append(genetuple)
    return genelist 

def oppredparser(opdoc):
    opdata = []
    with open(opdoc, "r") as file:
        for line in file.readlines():
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")
            currentline = tuple(currentline.split(","))
            opdata.append(currentline)
    op_dict = {gene: opn for opn, gene in opdata}
    return op_dict

def gbkreader(gbkdoc):
    convertdict = {"aaur":"AauR","algr":"AlgR","ampr":"AmpR","amrz":"AmrZ","anr":"Anr","argr":"ArgR","atzr":"AtzR","brlr":"BrlR","catr":"CatR","cifr":"CifR","colr":"ColR","cra":"Cra","cuer":"CueR","dest":"DesT","dnr":"DNR","exsa":"ExsA","fleq":"FleQ","fur":"Fur","hexr":"HexR","hrpl":"HrpL","ihf":"IHF","iscr":"IscR","lasr":"LasR","lexa":"LexA","mext":"MexT","mvfr":"MvfR","narl":"NarL","ntrc":"NtrC","ohrr":"OhrR","oxyr":"OxyR","pa2206":"PA2206","pcar":"PcaR","pchr":"PchR","phhr":"PhhR","phob":"PhoB","phzr":"PhzR","pltr":"PltR","psra":"PsrA","ptxr":"PtxR","ptxs":"PtxS","pvds":"PvdS","pyer":"PyeR","rcsb":"RcsB","rhpr":"RhpR","rpon":"RpoN","rsal":"RsaL","todt":"TodT","trpi":"TrpI","vfr":"Vfr","vqsm":"VqsM","vqsr":"VqsR","zur":"Zur"}
    with open(gbkdoc, "r") as file:
        namedgenes = [] #starts nested list for operon constituent genes in format [[LOCUSTAG1, OPERON1],...
        genefound = False
        for line in file:
            if "/gene=" in line:
                genereader = re.findall("/gene=\"(.*)\"", line)[0]
                if genereader in convertdict.keys():
                    geneadd = convertdict[genereader]
                    genefound = True
            if "locus_tag" in line and genefound == True:
                tagreader = re.findall("locus_tag=\"(.*)\"", line)[0]
                tmp = (geneadd, tagreader)
                namedgenes.append(tmp)
                genefound = False
    return namedgenes

def tf_igr_associator(genomeigrlist, currentgenome, geneigrassoc):
    tf_IGR_dict = defaultdict(list)
    for genome,locus,tf in genomeigrlist:
        if genome == currentgenome:
            try:
                tf_IGR_dict[tf].append(geneigrassoc[locus]) #end result should be {tf1:[igr1,igr2...]tf2:[igr5]...}
            except KeyError:
                if locus != "N/A-IGRin/isdoubleterminator":
                    print(f"ERROR - An IGR from table file for {tf} wasn't found in IGR_Genome_Gene file")
    return tf_IGR_dict
                
            
##########################
##########################
##########BODY############
##########################
##########################
datasetprompt = input("Which dataset are you updating the graphs for?('partial'/'comp'/'both' for datasets based upon partial/complete/complete+partial genomes respectively):")
if datasetprompt == "both":
    operonpredloc = "Data\\Processed_data\\operon_predictions_both"
    prokkaloc = "Data\\Processed_data\\prokka_out_both_genomes"
    graphfolder = "Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_bothgenomes"
    outputgraphfolder = "Results\\Graphml_full_regulatory_network_bothgenomes\\Regulatory_network_graph_"
if datasetprompt == "comp":
    operonpredloc = "Data\\Processed_data\\operon_predictions_comp"
    prokkaloc = "Data\\Processed_data\\prokka_out_complete_genomes"
    graphfolder = "Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_compgenomes"
    outputgraphfolder = "Results\\Graphml_full_regulatory_network_compgenomes\\Regulatory_network_graph_"
if datasetprompt == "partial":
    operonpredloc = "Data\\Processed_data\\operon_predictions_partial"
    prokkaloc = "Data\\Processed_data\\prokka_out_partial_genomes"
    graphfolder = "Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_partialgenomes"
    outputgraphfolder = "Results\\Graphml_full_regulatory_network_partialgenomes\\Regulatory_network_graph_"
igrgeneloc = "Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt"
hmmertabloc = "Data\\Hmmer_data\\nhmmerout"
operonpredictions = os.listdir(operonpredloc)
hmmerpreds = os.listdir(hmmertabloc)
prokkafiles = os.listdir(prokkaloc)


linedata = []
with open(igrgeneloc, "r") as file:
    for line in file.readlines():
        line = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","") #breaks down line into csv - IGR,GENOME,LOCUSTAG
        line = tuple(line.split(",")) #breaks line up by comma and groups values into tuple (IGR,GENOME,LOCUSTAG)
        if len(line) == 3: #Ensures that each line has 3 variables - IGR, genome and locustag, if not then there was an error in the IGR-GENOME-GENE associator
            linedata.append(line)
        else:
            print(f"Potential data issue, {line=}") #Reports error and line it occurred on, if it occurred
igr_dict = {gene:igr for igr,genome,gene in linedata} #forms a dictionary where genes are keys to the IGR theyre associated with


loopcounter = 0 #simple counter so progress in following loop can be reported in real time
genomeigrneighbour = tabreader(hmmerpreds, hmmertabloc) #returns list of tuples [(genome,locustag,transcriptionfactor),...], each tuple contains information to find IGR clusterID and the TF that binds that IGR
for x in prokkafiles: #loops for every folder from the prokka processed genomes
    loopcounter = loopcounter + 1
    tf_igr_interactions = tf_igr_associator(genomeigrneighbour, x, igr_dict) #returns the IGRs in genome 'x' that tfs bind to in form tf:igr1,igr2...
    regulatoryoperons = defaultdict(list) #Stores operons in genome 'x' that express regulators in format {tf1:[operon1,operon2...]}
    currentprediction = str(operonpredloc + "\\Operon_prediction_" + x + ".txt") #string concatenation to provide location of the operon prediction file for genome 'x'
    currentgbk = str(prokkaloc + "\\" + x + "\\" + x + ".gbk") #string concatenation to provide location of .gbk file for genome 'x'
    currentgmlloc = str(graphfolder + "\\IGR_OPERON_GML-" + x + ".graphml") #string concatenation to provide location of graphml file for genome 'x'
    gene_operon_dict = oppredparser(currentprediction) #returns dict of operons keyed by genes belonging to them - {gene1:op1, gene2:op1, gene3:op2...}
    operon_TFs = gbkreader(currentgbk) #returns list of tuples of transcription factors in genome and locus tag encoding it in format[(name,gene1),(name,gene2)...]
    for tf,tag in operon_TFs: #loops for every named gene in the gbk file
        try:
            regulatoryoperons[tf].append(gene_operon_dict[tag]) #finds operon the locus_tag(a.k.a gene) encoding a TF belongs to, makes dictionary of said operons keyed by TF they contain
        except KeyError: #Incase the gbkreader somehow returned a locustag that doesnt exist in genome 'x'
            if tag not in gene_operon_dict.keys():
                print(f"ERROR - {tag} is not a locus tag found in {x}")
    originalgraph = nx.read_graphml(currentgmlloc) #gets the data from the graphml files of genome 'x' which only has edges *from* IGR *to* the operon(s) it regulates
    for key in regulatoryoperons.keys(): #loops for each transcription factor known to be encoded in genome 'x'
        for operon in regulatoryoperons[key]: #loops for each operon associated with the TF(key) the above loop is currently on
            for cluster in tf_igr_interactions[key]: #loops for every IGR known to be associated with this TF in genome 'x'
                originalgraph.add_edge(operon, cluster) #for the current operon as defined by loop, and current IGR as defined by loop, forms an edge from operon to IGR
    nx.write_graphml(originalgraph, str(outputgraphfolder+x+".graphml")) #Writes a new graph in graphml format to a different folder than original
    print(f"The task is {(loopcounter/len(prokkafiles))*100} percent complete") 