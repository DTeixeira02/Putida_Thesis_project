#GRN_Statistics_calculator
#Determines average percent reduction in GRN node count following GRN processing by Motif_comparator_V2.4.py
#Also returns IGR counts, Operon counts and edge counts before and after merging for each graph

import networkx as nx
import os
import re
from collections import defaultdict
import statistics

versionprompt = 0
validchoice = False
validversion = False
while validchoice == False:
    datasetprompt = input("What dataset are you wanting to group IGRs for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH" or datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        validchoice = True
while validversion == False:
    versionprompt = input("Are you working with V2 or V3 complete graphs?('2' or '3') - V2 assumes only the hmmer specified occurrences of IGR are TF binding, V3 assumes hmmer specified occurrence of IGR AND all other occurrences of IGR is TF-binding site:")
    if versionprompt == "2" or versionprompt == "3":
        validversion = True


if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
    origgraphloc = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_Graphml_full_regulatory_network_bothgenomes"
    mastgraphloc = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_both_mastupdated"
    graphstatsout = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_both_stats.txt"
elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
    origgraphloc = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_Graphml_full_regulatory_network_compgenomes"
    mastgraphloc = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_comp_mastupdated"
    graphstatsout = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_comp_stats.txt"

ufids = os.listdir(mastgraphloc)
ids = []
for ufid in ufids:
    id = re.findall("GRN_MASTUPDATED_(.*).graphml",ufid)[0]
    ids.append(id)

before_igr_op_edge_counts = defaultdict(list)
after_igr_op_edge_counts = defaultdict(list)
befigrcounts = []
aftigrcounts = []
befopcounts = []
aftopcounts = []
befedgecounts = []
aftedgecounts = []
loopcounter = 0
for id in ids:
    loopcounter = loopcounter + 1
    print(loopcounter)
    before_igr_count = 0
    after_igr_count = 0
    before_operon_count = 0
    after_operon_count = 0
    before_edge_count = 0
    after_edge_count = 0
    originalgraph = nx.read_graphml(str(origgraphloc+"\\Regulatory_network_graph_"+id+".graphml"),node_type=str)
    mastgraph = nx.read_graphml(str(mastgraphloc+"\\GRN_MASTUPDATED_"+id+".graphml"),node_type=str)
    for node in list(originalgraph.nodes()):
        if "Cluster_" in node:
            before_igr_count = before_igr_count + 1
        elif "Operon" in node:
            before_operon_count = before_operon_count + 1
    before_edge_count = len(list(originalgraph.edges()))
    for node in list(mastgraph.nodes()):
        if "Cluster_" in node:
            after_igr_count = after_igr_count + 1
        elif "Operon" in node:
            after_operon_count = after_operon_count + 1
    after_edge_count = len(list(mastgraph.edges()))
    before_igr_op_edge_counts[id].extend([before_igr_count,before_operon_count,before_edge_count])
    after_igr_op_edge_counts[id].extend([after_igr_count,after_operon_count,after_edge_count])
    befigrcounts.append(before_igr_count)
    befopcounts.append(before_operon_count)
    befedgecounts.append(before_edge_count)
    aftigrcounts.append(after_igr_count)
    aftopcounts.append(after_operon_count)
    aftedgecounts.append(after_edge_count)
percentreductions = [(((befigr-aftigr)/befigr)*100) for befigr,aftigr in zip(befigrcounts,aftigrcounts)]
loopcounter = 0
with open(graphstatsout,"w") as file:
    file.write("This file lists the average percent reduction in IGR count as a result of IGR merges based on MAST data\n")
    file.write("#####################################################################\n")
    file.write("AVERAGE STATS ACROSS GENOMES\n")
    file.write(f"Mean IGR count before MAST merges: {sum(befigrcounts)/len(befigrcounts)} +/- {statistics.stdev(befigrcounts,(sum(befigrcounts)/len(befigrcounts)))}\n")
    file.write(f"Mean IGR count after MAST merges: {sum(aftigrcounts)/len(aftigrcounts)} +/- {statistics.stdev(aftigrcounts,(sum(aftigrcounts)/len(aftigrcounts)))}\n")
    file.write(f"Mean percent IGR count reduction: {sum(percentreductions)/len(percentreductions)} +/- {statistics.stdev(percentreductions,(sum(percentreductions)/len(percentreductions)))}\n")
    file.write(f"Mean Operon count: {sum(befopcounts)/len(befopcounts)} +/- {statistics.stdev(befopcounts,(sum(befopcounts)/len(befopcounts)))}\n")
    file.write(f"Mean Edge count: {sum(befedgecounts)/len(befedgecounts)} +/- {statistics.stdev(befedgecounts,(sum(befedgecounts)/len(befedgecounts)))}\n")
    file.write("#####################################################################\n")
    for id in before_igr_op_edge_counts.keys():
        file.write(f"{id}\n")
        file.write(f"IGR count before: {before_igr_op_edge_counts[id][0]}\n")
        file.write(f"IGR count after: {after_igr_op_edge_counts[id][0]}\n")
        file.write(f"Percent IGR node reduction: {percentreductions[loopcounter]}\n")
        file.write(f"Operon count: {before_igr_op_edge_counts[id][1]}\n")
        file.write(f"Edge count: {before_igr_op_edge_counts[id][2]}\n")
        file.write("#####################################################################\n")
        loopcounter = loopcounter + 1
