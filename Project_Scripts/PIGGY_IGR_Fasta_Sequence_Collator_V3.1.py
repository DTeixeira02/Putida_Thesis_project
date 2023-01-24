#PIGGY_IGR_Fasta_Sequence_Collator_V3.1(Output per genome)
#V3.1 differs from V3 in that V3.1 only uses neighbouring locus tags in fasta header, fixing length issue in MAST

import pandas as pd
from collections import defaultdict

def colreader(csvdata, colnumber):
    occurrencelist = []
    cellcounter = 0
    clusterIDs = csvdata.iloc[:,0]
    coldata = csvdata.iloc[:,colnumber]
    for cell in coldata:
        if pd.isnull(cell) == False:
            splittingdoubles = cell.split()
            for item in splittingdoubles:
                splitstring = item.split("_+_+_")
                header = str(splitstring[1]+"_+_+_"+splitstring[2]+"_+_+_"+splitstring[3])
                occurrencelist.append((clusterIDs[cellcounter],item,header))
        cellcounter = cellcounter + 1
    return occurrencelist

def IGRseqfilereader(sequences):
    occurrence_sequencedict = {}  
    with open(sequences, "r") as file2:
        startrecording = False
        for line in file2:
            writethisline = True
            formattedline = line.strip(">\n")
            if "GCA_" in formattedline:
                currentoccurrence = formattedline
                startrecording = True #meaning that the subsequent lines in the file will be the IGR sequence
                writethisline = False #prevents the current line(the header) from being included as part of the sequence - instead current line will be the key to the dict containing the sequence           
            if startrecording == True and writethisline == True:
                occurrence_sequencedict[currentoccurrence] = formattedline
    return occurrence_sequencedict


validchoice = False
while validchoice == False:
    datasetprompt = input("What piggy IGR dataset are you making GFF3 files for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        sequencesloc =  "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_sequences.txt"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\FASTA_IGRs_both_genomes_shortenedheaders"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_complete_genomes\\IGR_presence_absence.csv"
        sequencesloc =  "Internship_data\\Data\\Processed_data\\piggy_out_complete_genomes\\IGR_sequences.txt"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\FASTA_IGRs_comp_genomes_shortenedheaders"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")

presenceabsencecsv = pd.read_csv(presenceabsenceloc, low_memory=False)
IGRseqs = IGRseqfilereader(sequencesloc)
columniterator = 0
for col in presenceabsencecsv.columns:
    if "GCA_" in col:
        print(f"Currently collating sequences for {col} - Total task completion:{(columniterator/len(presenceabsencecsv.columns)*100)} %")
        IGRoccurrences = colreader(presenceabsencecsv, columniterator)
        with open(str(outputloc+"\\"+col+"_IGRs.fna"), "w") as file:
            for cluster,occurrence,header in IGRoccurrences:
                file.writelines(f">{header}\t{cluster}\n")
                file.writelines(f"{IGRseqs[occurrence]}\n")
    columniterator = columniterator + 1