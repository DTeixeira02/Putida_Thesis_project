#Motif_comparator_V2.3
#Motif_Comparator
#This is powered by spite, sleep-deprivation and black magic, dont ask me how it works
#If you run into issues, perform a ritual sacrifice to your processor and try again
#V2.1 uses the remcorr mastout file, which ignores the motifs with an E-value >0.05
#V2.2 same as 2.1 but uses different interspersing region scoring system: considers total length disparity and number of regions, doesnt consider order
#V2.2 adds to motif comparison function - weighting of length/sharedmotifs/motiforder/interspersingregions all appear as variables at start of function used during calculations, eases tweaking to scoring values
#V2.3 same as 2.1 and 2.2 but changes the merge function - IGRs deemed similar by mutual partner only merged if IGR similar to >50% of the current IGRs to be merged

#            for node in secondordernodes: #For every subIGR identified, i.e KEY <-> IGR <-> subIGR <> subsubigr
#                sharedcount = 0 #Counts the number of subsubIGRs found 
#                for subigr in igrcomparisondict[node]: #for every subsubigr found similar to a subigr found similar to the IGR found directly similar to the initial key
#                    if subigr in nodelist:
#                        sharedcount = sharedcount + 1
#                if (sharedcount/len(igrcomparisondict[node])) >= 0.7: #if over 70% of the IGRs found under subsubIGR
#                    for subigr in igrcomparisondict[node]:
#                        if subigr not in alreadymergedlist and subigr not in nodelist:
#                            nodelist.append(subigr)
#                if node in alreadymergedlist:
#                    print(f"####   {value} is already merged with a different cluster than {key}   ####")


import os
from collections import defaultdict
import re
import networkx as nx

def merge_nodes(G,orignodes):
    attemptmerge = True
    nodes = []
    for node in orignodes:
        if node in list(G.nodes()):
            nodes.append(node)
    new_node = ','.join(nodes)
    if len(nodes) <= 1:
        attemptmerge = False
    if attemptmerge == True:
        G.add_node(new_node) # Add the 'merged' node
        edgelist = list(G.edges())
        for n1,n2 in edgelist:
            if n1 in nodes:
                G.add_edge(new_node,n2)
            elif n2 in nodes:
                G.add_edge(n1,new_node)
        for n in nodes: # remove the merged nodes
            G.remove_node(n)
        #print(f"Merged node {new_node} created from {orignodes} - Omitted nodes were double terminators")
    return G,nodes

def mastreader(masttxt):
    occurrencemotifs = [] #nested list - each contained list in format [Tag1_+_+_tag2_, clusterID, Motifs]
    locustoclusterdict = {}
    motifdict = defaultdict(list)
    cluster_to_diagrams_dict = defaultdict(list)
    datamotiffound = False #Following False variables used to identify which section of MAST txt currently being read.
    scoringseqfound = False
    motifdiagramfound = False
    significant = False #Used during motifdiagramfound section to identify if Unique IGR occurrence ID has been found - i.e the tag_+_+_tag_+_+_orientation string
    with open(masttxt, "r") as file:
        for line in file:
            if "DATABASE AND MOTIFS" in line:
                datamotiffound = True
            if "PAIRWISE MOTIF CORRELATIONS:" in line:
                datamotiffound = False
            if "SECTION I: HIGH-SCORING SEQUENCES" in line:
                scoringseqfound = True
            if "SECTION II: MOTIF DIAGRAMS" in line:
                scoringseqfound = False
                motifdiagramfound = True
            if datamotiffound == True:
                if "MEME-" in line:
                    splitline = line.split()
                    if "MEME-" in splitline[2]:
                        motifdict[splitline[0]] = [splitline[1],splitline[2],splitline[3],splitline[4]] #Dictionary with motif number as key, list in format motifID, ALT-ID, Width, Best possible match/consensus sequence?
            if scoringseqfound == True:
                if "_+_+_" in line:
                    splitline = line.split()
                    locustag = splitline[0].split("_+_+_")[0]
                    locustoclusterdict[locustag] = splitline[1]
            if motifdiagramfound == True:
                if "_+_+_" in line:
                    splitline = line.split()
                    evalpower = splitline[1].split("e")
                    if len(evalpower) == 2:
                        if int(evalpower[1]) > -5: #checks that evalue is -5 or less(i.e more negative) as i thought it would be an ok cutoff
                            significant = False
                        elif int(evalpower[1]) <= -5:
                            significant = True
                    if significant == True and len(evalpower) == 2:
                        locustag = splitline[0].split("_+_+_")[0]
                        try:
                            occurrencemotifs.append([splitline[0],locustoclusterdict[locustag],splitline[2]])
                            cluster_to_diagrams_dict[locustoclusterdict[locustag]].append(splitline[2])
                        except KeyError:
                            print(f"No cluster associated with {locustag}")
            if "SECTION III: ANNOTATED SEQUENCES" in line:
                break
    return motifdict,occurrencemotifs,cluster_to_diagrams_dict

def motifcomparison(IGRoccs): #takes all occurrences of IGRs and their 'motif structure' and compares them to decide which structures are similar enough to cluster, and therefore which IGRs are similar enough to cluster
    similaroccurrences = [] #Nested list - format [CurrentoccCluster,OtheroccCluster,CurrentoccTags,OtheroccTags]
    intersperseweighting = 125
    sharemotifweighting = 300
    shareorderweighting = 300
    lengthdispweighting = 250
    totalthresh = 725
    minthresh = 650
    for occurrence in IGRoccs:
        currentoccmotifs = []
        currentoccinterspersing = []
        currentoccurrence = occurrence
        currentoccdiagram = re.findall("\d+|\[\+\d?\d?\]|\[\-\d?\d?\]", currentoccurrence[2]) #this abomination checks for the numbers inside the brackets- \d means any number 0-9, ? means 0 or 1 occurrence, | means either or(used since sign can be pos or neg). Going per section separated by |, sec1 checks for the value in the very 1st pos assuming its not motif, second section does the same but for very last
        for segment in currentoccdiagram:
            if "[" in segment:
                currentoccmotifs.append(segment)
            if "[" not in segment:
                currentoccinterspersing.append(segment)
        currentoccinterlength = sum(map(int,currentoccinterspersing))
        for otherocc in IGRoccs:
            attemptcomparison = True
            if currentoccurrence[1] == otherocc[1]:
                #print(f"Not attempting to compare Current Occurrence {currentoccurrence[0]} to Other occurrence {otherocc[0]} as they are occurrences of the same IGR {currentoccurrence[1]}")
                attemptcomparison = False
            if attemptcomparison == True:
                sharedcounter = 0
                otheroccmotifs = [] #contains only the motifs as they appear in otherocc diagram - NO INVERTED SIGN MOTIFS
                otheroccinterspersing = [] #contains only the interspersing regions as they appear in otherocc diagram
                otheroccinterspersingdupe = []
                otheroccmotifsandinverse = [] #Contains all motifs as they appear in otherocc diagram, aswell as the same motifs with the inverted sign  #Due to how the code works, motifs will be removed from other occ motifs, but they will still be needed later, hence 2 copies
                sharedmotifscurrentordered = []
                sharedmotifsotherordered = []
                sharedmotifsandinverse = []
                otheroccdiagram = re.findall("\d+|\[\+\d?\d?\]|\[\-\d?\d?\]", otherocc[2])
                
                for segment in otheroccdiagram:
                    if "[" in segment:
                        otheroccmotifs.append(segment)
                        if "+" in segment:
                            oppstrand = segment.replace("+","-")
                            otheroccmotifsandinverse.extend([segment,oppstrand]) 
                        elif "-" in segment:
                            oppstrand = segment.replace("-","+")
                            otheroccmotifsandinverse.extend([segment,oppstrand])
                        else:
                            print(f"Error in {segment}")
                    if "[" not in segment:
                        otheroccinterspersing.append(segment)
                        otheroccinterspersingdupe.append(segment)
               
               
                if len(currentoccmotifs) > len(otheroccmotifs):
                    lengthcontrib = lengthdispweighting*(len(otheroccmotifs)/len(currentoccmotifs))
                if len(currentoccmotifs) <= len(otheroccmotifs): #The equal sign included here but not in the other to allow cases where len is equal, but to prevent it doing it twice - wouldnt cause error, just for performance
                    lengthcontrib = lengthdispweighting*(len(currentoccmotifs)/len(otheroccmotifs))
               
               
                otheroccinterlength = sum(map(int,otheroccinterspersing))
                if len(currentoccinterspersing) != 0 and len(otheroccinterspersing) != 0:
                    if len(currentoccinterspersing) >= len(otheroccinterspersing):
                        if currentoccinterlength >= otheroccinterlength:
                            interspersecontrib = intersperseweighting*((otheroccinterlength/currentoccinterlength)*(len(otheroccinterspersing)/len(currentoccinterspersing)))
                        if currentoccinterlength < otheroccinterlength:
                            interspersecontrib = intersperseweighting*((currentoccinterlength/otheroccinterlength)*(len(otheroccinterspersing)/len(currentoccinterspersing)))
                    elif len(currentoccinterspersing) < len(otheroccinterspersing):
                        if currentoccinterlength >= otheroccinterlength:
                            interspersecontrib = intersperseweighting*((otheroccinterlength/currentoccinterlength)*(len(currentoccinterspersing)/len(otheroccinterspersing)))
                        if currentoccinterlength < otheroccinterlength:
                            interspersecontrib = intersperseweighting*((currentoccinterlength/otheroccinterlength)*(len(currentoccinterspersing)/len(otheroccinterspersing)))
                if len(currentoccinterspersing) == 0 or len(otheroccinterspersing) == 0:
                    interspersecontrib = 0

  
                for motif in currentoccmotifs: #currentoccmotifs is only the motifs as they appear in currentoccdiagram
                    if motif in otheroccmotifsandinverse:
                        sharedcounter = sharedcounter + 1
                        if "+" in motif:
                            reversemotif = motif.replace("+","-")
                        elif "-" in motif:
                            reversemotif = motif.replace("-","+")
                        sharedmotifscurrentordered.append(motif)
                        sharedmotifsandinverse.extend([motif,reversemotif])
                        otheroccmotifsandinverse.remove(motif) #If the same motif occurs twice in currentocc but only once in otherocc, this removal prevents double counting of sharing that motif
                        otheroccmotifsandinverse.remove(reversemotif) #Same reason above, just removes partner motif with opposite sign              
                for motif in otheroccmotifs: #otheroccmotifs is only the motifs as they appear in otheroccdiagram
                    if motif in sharedmotifsandinverse: #sharedmotifsandinverse contains only the motifs that occurred in both alongside the same motif value w/opp sign
                        sharedmotifsotherordered.append(motif) #Appends the shared motifs in the order they occur in otherocc
                        if "+" in motif:
                            reversemotif = motif.replace("+","-")
                        elif "-" in motif:
                            reversemotif = motif.replace("-","+")
                        sharedmotifsandinverse.remove(motif)
                        sharedmotifsandinverse.remove(reversemotif)
                
                orderedcounterforward = sharedcounter #Counter deducted for each motif that isnt in the same order in otherocc vs currentocc
                orderedcounterreverse = sharedcounter #same as above, but used when currentocc compared to reverse of otherocc
                for i in range(len(sharedmotifscurrentordered)):
                    if "+" in sharedmotifscurrentordered[i]:
                        revcurrentmotif = sharedmotifscurrentordered[i].replace("+","-")
                    if "-" in sharedmotifscurrentordered[i]:
                        revcurrentmotif = sharedmotifscurrentordered[i].replace("-","+")
                    if sharedmotifscurrentordered[i] != sharedmotifsotherordered[i] and revcurrentmotif != sharedmotifsotherordered[i]:
                        orderedcounterforward = orderedcounterforward-1
                    if sharedmotifscurrentordered[i] != sharedmotifsotherordered[len(sharedmotifsotherordered)-1-i] and revcurrentmotif != sharedmotifsotherordered[len(sharedmotifsotherordered)-1-i]: #Compares 'I'th motif from the left of currentocc to the 'I'th motif from the right of otherocc
                        orderedcounterreverse = orderedcounterreverse-1

                if sharedcounter >2: #If shared counter is 1 or 2, the ordercontrib is max by default as it is always in ordeer with itself. 0 causes division errors    
                    if orderedcounterforward > orderedcounterreverse: #If statement to ensure the otherocc orientation with greatest order similarity to currentocc is selected
                        ordercontrib = shareorderweighting*(orderedcounterforward/sharedcounter) #Deduction based on proportion of unordered shared motifs to shared motifs total                  
                    if orderedcounterforward <= orderedcounterreverse:
                        ordercontrib = shareorderweighting*(orderedcounterreverse/sharedcounter)
                if sharedcounter == 0 or sharedcounter == 1 or sharedcounter == 2:
                    if (sharedcounter/len(currentoccmotifs)) >= 0.5 and (sharedcounter/len(otheroccmotifs)) >= 0.5: #If the shared motifs make up over half of the motifs in both IGRs
                        ordercontrib = 165 #Smaller amount due to fact sharing less motifs naturally means more are likely in order - lower order score makes the other factors be more important in final scoring
                    else:
                        ordercontrib = 0
                if len(currentoccmotifs) > 0 and len(otheroccmotifs) > 0:
                    sharecontrib = sharemotifweighting*((sharedcounter/len(currentoccmotifs))*(sharedcounter/len(otheroccmotifs)))
                if len(currentoccmotifs) == 0 or len(otheroccmotifs) == 0:
                    sharecontrib = 0
                totalscore = lengthcontrib+sharecontrib+ordercontrib+interspersecontrib
                if totalscore >= totalthresh:
                    similaroccurrences.append([currentoccurrence[1],otherocc[1],currentoccurrence[0],otherocc[0],currentoccurrence[2],otherocc[2],totalscore])
                #if totalscore >= minthresh and totalscore < totalthresh:
                #    print(f"{currentoccurrence[1]},{otherocc[1]}\t{currentoccurrence[2]}\t{otherocc[2]}\t{totalscore}")
    return similaroccurrences


print("This script requires mast.txt file for each genome you want to group IGRs for and the \"complete\" network graphs produced as detailed in program execution instructions. Please run mast first if not yet done")
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
manualreviewprompt = input("Do you wish to manually review IGR similarity for proposed similar IGRs?(Y/N):")

if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
    origgraphloc = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_Graphml_full_regulatory_network_bothgenomes"
    newgraphout = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_both_mastupdated"
    mergedclustersout = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_both_mastmerged_IGRs"
    mastloc = "Thesis_project_data\\Data\\Processed_data\\mastout_remcorr"
    if os.path.isdir(newgraphout) == False:
        os.mkdir(newgraphout)
    if os.path.isdir(mergedclustersout) == False:
        os.mkdir(mergedclustersout)
elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
    origgraphloc = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_Graphml_full_regulatory_network_compgenomes"
    newgraphout = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_comp_mastupdated"
    mergedclustersout = f"Thesis_project_data\\Data\\Processed_data\\V{versionprompt}_GRN_comp_mastmerged_IGRs"
    mastloc =  "Thesis_project_data\\Data\\Processed_data\\mastout_remcorr"
    if os.path.isdir(newgraphout) == False:
        os.mkdir(newgraphout)
    if os.path.isdir(mergedclustersout) == False:
        os.mkdir(mergedclustersout)
genomeids = os.listdir(mastloc)

for id in genomeids:    
    reviewed = False
    mastfile = str(mastloc+"\\"+id+"\\mast.txt")
    motifdict,occurrencemotifs,cluster_diagram_dict = mastreader(mastfile)
    print("###############################")
    print(f"Beginning similarity comparison for genome {id}")
    igrstomerge = motifcomparison(occurrencemotifs) #Contained lists in structure [CurrentoccCluster,OtheroccCluster,CurrentoccTags,OtheroccTags,motifdiag1,motifdaig2,score]
    if manualreviewprompt == "Y" or manualreviewprompt == "y": #CODE IS BUGGED HERE, DO NOT USE
        while reviewed == False:
            counter = 0
            for item in igrstomerge:
                counter = counter + 1
                if item[6] <= 740:
                    print(f"{counter}. {item[0]} {item[6]}   vs   {item[1]}\n{item[4]}\n{item[5]}\n")
            print("###############################")
            print("Input the number associated to the comparisons you believe aren't similar, separated by a whitespace - They will be deleted")
            print("If you don't want to delete any comparisons, press enter without entering a value")
            deleteprompt = input("Comparisons to delete:")
            itemstodelete = deleteprompt.split()
            for deletion in itemstodelete:
                igrstomerge.pop((deletion-1))
                print(f"Deleted {deletion}. {igrstomerge[deletion-1][0]} {igrstomerge[deletion-1][4]} vs {igrstomerge[deletion-1][1]} {igrstomerge[deletion-1][5]} Total score:{igrstomerge[deletion-1][6]}")
            endprompt = input("Are you finished deletions?(Y/N):")
            if endprompt == "Y" or endprompt == "y":
                reviewed = True
    if manualreviewprompt == "N" or manualreviewprompt == "n":
        print("User selected not to manually review IGR similarity")

    originalgraph = nx.read_graphml(str(origgraphloc+"\\Regulatory_network_graph_"+id+".graphml"),node_type=str)
    igrcomparisondict = defaultdict(list)
    for entry in igrstomerge:
        if entry[0] in igrcomparisondict.keys():
            if entry[1] not in igrcomparisondict[entry[0]]:
                igrcomparisondict[entry[0]].append(entry[1])
        if entry[1] in igrcomparisondict.keys():
            if entry[0] not in igrcomparisondict[entry[1]]:
                igrcomparisondict[entry[1]].append(entry[0])
        if entry[0] not in igrcomparisondict.keys():
            igrcomparisondict[entry[0]].append(entry[1])
        if entry[1] not in igrcomparisondict.keys():
            igrcomparisondict[entry[1]].append(entry[0])
    nodelist = []
    nodestomerge = [] #nested list storing each cluster that should be merged
    alreadymergedlist = []
    print("Beginning cluster merges")
    for key in igrcomparisondict.keys():
        if key not in alreadymergedlist:
            nodelist = [key]
            for value in igrcomparisondict[key]:
                if value not in alreadymergedlist and value not in nodelist:
                    nodelist.append(value)
            firstordernodes = nodelist
            secondordernodes = []
            for node in firstordernodes: #For every IGR similar to the original Key
                if node not in alreadymergedlist and node != key: #If those IGRs havent been merged already and arent the original key
                    if node in igrcomparisondict.keys(): #If those IGRs have a similarity to atleast one other IGR(Should be guaranteed)
                        for node2 in igrcomparisondict[node]: #For every subIGR similar to those IGRs found similar to the key, note the subIGRs arent necessarily deemed similar to the key
                            if node2 not in alreadymergedlist and node2 not in nodelist and node2 != key: #If the subIGR isnt already merged, isnt the key, and isnt already part of nodelist
                                nodelist.append(node2) 
                                secondordernodes.append(node2)
                    if node not in igrcomparisondict.keys():
                        print(f"Potential error - {node} somehow doesn't have any similarites despite appearing under {key}")
            alreadymergedlist.extend(nodelist)
            if len(nodelist) > 1:
                nodestomerge.append(nodelist)
            elif len(nodelist) <= 1:
                print(f"All clusters similar to {key} have already been merged with other clusters, {key} will not be merged.")
    mergenodelist=[]
    for nodelist in nodestomerge:
        originalgraph,mergednodes = merge_nodes(originalgraph,nodelist)
        if len(mergednodes) > 1:
            mergenodelist.append(mergednodes)
    print("Merges complete - Writing graph to file")
    nx.write_graphml(originalgraph, str(newgraphout+"\\GRN_MASTUPDATED_"+id+".graphml"))
    print("Graphs written - writing merged IGR list to text file")
    with open(str(mergedclustersout+"\\Merged_IGRs_"+id+".txt"), "w") as file:
        dupe_checking = []
        for merge in mergenodelist:
            clusters_with_diagram = []
            for igr in merge:
                if igr in dupe_checking:
                    print(f"{igr} is merged more than once in {id}")
                if len(cluster_diagram_dict[igr]) > 1:
                    diagram = "|".join(cluster_diagram_dict[igr])
                    formattedentry = str(igr+"|"+diagram)
                    clusters_with_diagram.append(formattedentry)
                if len(cluster_diagram_dict[igr]) <= 1:
                    formattedentry = str(igr+"|"+cluster_diagram_dict[igr][0])
                    clusters_with_diagram.append(formattedentry)
                dupe_checking.append(igr)
            line = '\t'.join(clusters_with_diagram)
            file.write(f"{line}\n\n")