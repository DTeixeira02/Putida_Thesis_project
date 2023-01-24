#IGR_gene_associator_V2
import os
import pandas as pd
import re
from collections import defaultdict

def blankdel(cell): #function for removing empty cells - this doesnt delete them from the base file, just the rowdata list, errors occurred when blanks present: easier to fix by removing blanks
    if cell == (""):
        return False
    return True  

def presenceabsencecsvreader(csvdata):
    IGR_to_tag_dict = defaultdict(list)
    doubleterminators = defaultdict(list)
    multiassocs = []
    cellcounter = 0
    rowcounter = 0
    clusterIDs = csvdata.iloc[:,0]
    denominator = len(clusterIDs)
    for cluster in clusterIDs:
        rowdata = csvdata.iloc[rowcounter,:]
        for cell in rowdata:
            if pd.isnull(cell) == False:
                if "GCA_" in str(cell):
                    splittingdoubles = cell.split()
                    for item in splittingdoubles:
                        doubleassoc = False
                        splitstring = item.split("_+_+_")
                        if splitstring[3] != "DT" and (splitstring[1] in IGR_to_tag_dict.keys() or splitstring[2] in IGR_to_tag_dict.keys()):
                            doubleassoc = True
                            print(f"Potential error, {splitstring[1]} or {splitstring[2]} associated with multiple clusters")
                            if splitstring[3] == "CO_F":
                                multiassocs.append([cluster,splitstring[2],splitstring[0]])
                            if splitstring[3] == "CO_R":
                                multiassocs.append([cluster,splitstring[1],splitstring[0]])
                            if splitstring[3] == "DP":
                                multiassocs.append([cluster,splitstring[2],splitstring[0]])
                                multiassocs.append([cluster,splitstring[1],splitstring[0]])
                        if splitstring[3] == "CO_F" and doubleassoc == False:
                            IGR_to_tag_dict[cluster].append([splitstring[2],splitstring[0]])
                        if splitstring[3] == "CO_R" and doubleassoc == False:
                            IGR_to_tag_dict[cluster].append([splitstring[1],splitstring[0]])
                        if splitstring[3] == "DP" and doubleassoc == False:
                            IGR_to_tag_dict[cluster].append([splitstring[2],splitstring[0]])
                            IGR_to_tag_dict[cluster].append([splitstring[1],splitstring[0]])
                        if splitstring[3] == "DT" and doubleassoc == False:
                            doubleterminators[cluster].append([splitstring[1],splitstring[2],splitstring[0]]) 
            cellcounter = cellcounter + 1
        rowcounter = rowcounter+1
        print(f"Task at {(rowcounter/denominator)*100} \% completion")
    return IGR_to_tag_dict,doubleterminators,multiassocs


validchoice = False
while validchoice == False:
    datasetprompt = input("What piggy IGR dataset are you making IGR-Gene association files for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        outputloc = "Thesis_project_data\\Data\\Irrelevant_data\\IGR_gene_associations\\IGR_assocs_bothgenomes.txt"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        outputloc = "Thesis_project_data\\Data\\Irrelevant_data\\IGR_gene_associations\\IGR_assocs_compgenomes.txt"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")


presenceabsencecsv = pd.read_csv(presenceabsenceloc, low_memory=False)
tagigrdict, dtdict, multiassoclist = presenceabsencecsvreader(presenceabsencecsv)
denominator = len(tagigrdict.keys())
loopcounter = 0
print("###############################################################")
print("Beginning writing process")
print(f"Task at {(loopcounter/denominator)*100}\%")
with open(outputloc, "w") as file:
    for cluster in tagigrdict.keys():
        for gene,genome in tagigrdict[cluster]:
            file.write(f"{cluster}\t{gene}\t{genome}\n")
        for igr,gene,genome in multiassoclist:
            if igr == cluster:
                file.write(f"{cluster}\t{gene}\t{genome}\n")
        for gene1,gene2,genome in dtdict[cluster]:
            file.write(f"{cluster}\tDOUBLE_TERMINATOR_BETWEEN_{gene1}_{gene2}\t{genome}\n")
        print(f"Task at {(loopcounter/denominator)*100}\% completion")
        loopcounter = loopcounter+1