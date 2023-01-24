import pandas as pd
import os
import re
from collections import defaultdict

#gff3 format: ID_of_original_genomefragment/wholegenome  PIGGY    Intergenic Region    start    end    .    .    .    .

def stringreader(currentcolumnnumber, dataframe, taglocations):
    features = []
    rowiterator = 0
    clusterIDs = dataframe.iloc[:,0]
    columndata = dataframe.iloc[:,currentcolumnnumber]
    for cell in columndata:
        if pd.isnull(cell) == False:
            splittingdoubles = cell.split()
            if len(splittingdoubles) > 1:
                for occurrence in splittingdoubles:
                    dividedstring = occurrence.split("_+_+_")
                    features.append([taglocations[dividedstring[1]][0],"PIGGY","IGR",int(taglocations[dividedstring[1]][2])+1,int(taglocations[dividedstring[2]][1])-1,".",".",".",str("ID="+clusterIDs[rowiterator]+";upstream="+dividedstring[1]+";downstream="+dividedstring[2]+";direction="+dividedstring[3])])
                    if taglocations[dividedstring[1]][0] != taglocations[dividedstring[2]][0]:
                        print("IGR between genes on different genome fragments of partial genome")
            else:
                dividedstring = cell.split("_+_+_")
                features.append([taglocations[dividedstring[1]][0],"PIGGY","IGR",int(taglocations[dividedstring[1]][2])+1,int(taglocations[dividedstring[2]][1])-1,".",".",".",str("ID="+clusterIDs[rowiterator]+";upstream="+dividedstring[1]+";downstream="+dividedstring[2]+";direction="+dividedstring[3])])
        rowiterator = rowiterator + 1 #used to keep track of current row we're on so that the correct cluster is associated to the current IGR occurrence
    return features

def tagtolocation(currentgenome, genegff3loc):
    sequenceregions = []
    fastas = defaultdict(list)
    tmp = []
    tmp2 = []
    tagdict = defaultdict(list)
    fastafound = False
    validID = False
    with open(str(genegff3loc+"\\"+currentgenome+".gff"), "r") as file2:
        for line in file2.readlines():
            if "ID=" in line: #needed to make sure only feature lines are scanned - errors with refindall otherwise
                splitstring = line.split()
                attributestring = splitstring[8].split(";")
                locustag = re.findall("ID=(.*)", attributestring[0])[0] #scans feature string for locus tag
                tagdict[locustag] = [splitstring[0],splitstring[3],splitstring[4]] #format: SeqID, Start pos, end pos
                if splitstring[0] not in tmp:
                    tmp.append(splitstring[0])
            if "##sequence-region" in line:
                regionsplit = line.split()
                sequenceregions.append(regionsplit)
                tmp2.append(regionsplit[1])
            if "##FASTA" in line:
                fastafound = True
            if fastafound == True:
                writethisline = True
                scanstring = line.strip(">\n")
                if scanstring in tmp2:
                    currentID = scanstring
                    validID = True
                    writethisline = False 
                if validID == True and writethisline == True: #needed to prevent attempts at appending line before a currentID has been established, writethisline prevents the header being added to the list associated with the key
                    fastas[currentID].append(scanstring)
        for region in sequenceregions:
            if region[1] not in tmp:
                sequenceregions.remove(region)
    return tagdict, sequenceregions, fastas


##################################################################################################################################################################################

validchoice = False
while validchoice == False:
    datasetprompt = input("What piggy IGR dataset are you making GFF3 files for?('both' or 'comp'):")
    if datasetprompt == "both" or "Both" or "BOTH":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\GFF3_IGRs_both_genomes"
        genegff3loc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\gff_files"
        validchoice = True
    elif datasetprompt == "comp" or "Comp" or "COMP":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_complete_genomes\\IGR_presence_absence.csv"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\GFF3_IGRs_comp_genomes"
        genegff3loc = "Internship_data\\Data\\Processed_data\\piggy_out_complete_genomes\\gff_files"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset")


presenceabsencecsv = pd.read_csv(presenceabsenceloc, low_memory=False)
columniterator = 0
for col in presenceabsencecsv.columns:
    if "GCA_" in col:
        print("##########")
        print(f"GFF3 formatting IGR data for {col}")
        tagdict, seqregions, fastadict = tagtolocation(col, genegff3loc)
        gfffeatures = stringreader(columniterator, presenceabsencecsv, tagdict)
        with open(str(outputloc+"\\GFF3_"+col+".txt"), "w") as file:
            file.writelines("##gff-version 3\n")
            for region in seqregions:
                file.writelines(str(" ".join(region))+"\n")
            for feature in gfffeatures:
                file.writelines(f"{feature[0]}\t{feature[1]}\t{feature[2]}\t{feature[3]}\t{feature[4]}\t{feature[5]}\t{feature[6]}\t{feature[7]}\t{feature[8]}\n")
            file.writelines("##FASTA\n")
            for region in seqregions:
                file.writelines(f">{region[1]}\n")
                for nucleotideseq in fastadict[region[1]]:
                    file.writelines(f"{nucleotideseq}\n")
    columniterator = columniterator + 1
