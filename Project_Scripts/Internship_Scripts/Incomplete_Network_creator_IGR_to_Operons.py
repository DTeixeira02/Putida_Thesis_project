#Graphing_script-input:IGR_Operon_assoc
#THIS SCRIPT ONLY FORMS A GRAPH WITH EDGES *FROM IGR TO OPERONS*, THE OUTPUT MUST BE PROCESSED WITH ANOTHER SCRIPT TO OBTAIN FULL REGULATORY NETWORK

import networkx as nx
import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt

def parse_igr_pred(igrdoc, genid):
    igr_data = []
    with open(igrdoc, "r") as file:
        for line in file.readlines():
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")  
            currentline = tuple(currentline.split(","))  
            if currentline[1] == genid and len(currentline) == 3:  
                igr_data.append(currentline)
            elif len(currentline) != 3:  
                print(f"Potential data issue, {line=}")
    gene_igr_dict = {gene: igr for igr, genome, gene in igr_data}
    return gene_igr_dict
            
def parse_op_pred(opdoc):
    op_data = []
    operon_gene_dict = defaultdict(list)
    with open(opdoc, "r") as file1:
        for line in file1.readlines():
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")  
            currentline = tuple(currentline.split(","))
            if len(currentline) == 2:
                op_data.append(currentline)
            else:
                print(f"Potential data issue, {line=}")
    for opn, gene in op_data:
        operon_gene_dict[opn].append(gene)
    return operon_gene_dict

datasetprompt = input("Which dataset will you form incomplete graphs for?(partial/comp/both):")
if datasetprompt == "comp":
    igrloc = "Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt" #This is meant to say 'both_genomes' so that locus tags are consistent between datasets/graphs
    oploc = "Data\\Processed_data\\operon_predictions_comp"
    graphloc = "Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_compgenomes\\IGR_OPERON_GML-"
elif datasetprompt == "both":
    igrloc = "Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt"
    oploc = "Data\\Processed_data\\operon_predictions_both"
    graphloc = "Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_bothgenomes\\IGR_OPERON_GML-"
elif datasetprompt == "partial":
    igrloc = "Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt"
    oploc = "Data\\Processed_data\\operon_predictions_partial"
    graphloc = "Data\\Processed_data\\Incomplete_networks\\Graphml_incomplete_graphs_partialgenomes\\IGR_OPERON_GML-"

collection = os.listdir(oploc)

for i in range(len(collection)):
    operon_no_assocs = []
    op_igr = {}
    currentgenomeloc = str(oploc + "\\" + collection[i])
    currentgenomeid = str(re.findall("Operon_prediction_(.*).txt", collection[i])[0])
    geneigrdict = parse_igr_pred(igrloc, currentgenomeid) #returns genes associated with IGR in form {gene1: IGR1, gene2:IGR439824...}
    op_gene_dict = parse_op_pred(currentgenomeloc) #returns dictionary of gene lists keyed by the operon they belong to, i.e {opn1:[gene1,gene2,gene3],opn2:[gene4,gene5,gene6]}
    for operon in op_gene_dict.keys(): #repeats for every operon predicted for the genome
        if operon not in op_igr.keys(): #ensures that no operon will be associated with more than one cluster
            if op_gene_dict[operon][0] in geneigrdict.keys(): #checks if theres an IGR-gene association for the first gene in the operon
                op_igr[operon] = geneigrdict[op_gene_dict[operon][0]]
            elif op_gene_dict[operon][len(op_gene_dict[operon])-1] in geneigrdict.keys(): #if the above fails, checks if theres an IGR-gene association for the last gene in the operon
                op_igr[operon] = geneigrdict[op_gene_dict[operon][len(op_gene_dict[operon])-1]]
            else: #If no IGR was found for either of the above genes, then the operon number is recorded in a separate list so nodes are still included in the graph
                operon_no_assocs.append(operon)
    print("#############")
    print(f"Total no. operons = {len(op_gene_dict.keys())}")
    print(f"Total no. operons with IGRs = {len(op_igr.keys())}")
    print(f"Total no. operons without IGRs = {len(operon_no_assocs)}")
    if len(op_gene_dict.keys()) != (len(op_igr.keys())+len(operon_no_assocs)):
        print("ERROR - Some operons have either been skipped or double counted in the op_gene dict")
    igr_graph = nx.DiGraph()
    for operon, igr in op_igr.items():
        igr_graph.add_edge(igr, operon)
    for operon in operon_no_assocs:
        igr_graph.add_node(operon)
    nx.write_graphml(igr_graph, str(graphloc+currentgenomeid+".graphml"))  
    print(f"The task is {(i/len(collection))*100} percent done")