#PIGGY_IGR_Fasta_Sequence_Collator
import pandas as pd
from collections import defaultdict

def rowreader(rowdata, sequences):
    occurrence_sequencedict = defaultdict(list)
    occurrencelist = []
    for cell in rowdata:
        if pd.isnull(cell) == False:
            if "GCA_" in str(cell):
                splittingdoubles = cell.split()
                for item in splittingdoubles:
                    occurrencelist.append(item)
    with open(sequences, "r") as file2:
        startrecording = False
        for line in file2:
            writethisline = True
            formattedline = line.strip(">\n")
            if "GCA_" in formattedline:
                if formattedline in occurrencelist:
                    currentoccurrence = formattedline
                    startrecording = True #meaning that the subsequent lines in the file will be the IGR sequence
                    writethisline = False #prevents the current line(the header) from being included as part of the sequence - instead current line will be the key to the dict containing the sequence
                else:
                    startrecording = False            
            if startrecording == True and writethisline == True:
                occurrence_sequencedict[currentoccurrence].append(formattedline)
    return occurrence_sequencedict


validchoice = False
while validchoice == False:
    datasetprompt = input("What piggy IGR dataset are you making GFF3 files for?('both' or 'comp'):")
    if datasetprompt == "both" or "Both" or "BOTH":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
        sequencesloc =  "Internship_data\\Data\\Processed_data\\piggy_out_both_genomes\\IGR_sequences.txt"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\FASTA_IGRs_both_genomes"
        validchoice = True
    elif datasetprompt == "comp" or "Comp" or "COMP":
        presenceabsenceloc = "Internship_data\\Data\\Processed_data\\piggy_out_complete_genomes\\IGR_presence_absence.csv"
        sequencesloc =  "Internship_data\\Data\\Processed_data\\piggy_out_complete_genomes\\IGR_sequences.txt"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\FASTA_IGRs_comp_genomes"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset")

presenceabsencecsv = pd.read_csv(presenceabsenceloc, low_memory=False)
cluster_occurrence_dict = defaultdict(list)
for i in range(len(presenceabsencecsv.iloc[:,0])):
    currentrow = presenceabsencecsv.iloc[i,:]
    print(f"Currently collating sequences for {currentrow[0]} - Total task completion:{(i/len(presenceabsencecsv.iloc[:,0]))*100}\%")
    seqdict = rowreader(currentrow, sequencesloc)
    with open(str(outputloc+"\\"+currentrow[0]+"_collated_sequences.ffn"), "w") as file:
        for key in seqdict.keys():
            file.writelines(f">{key}\n")
            for seqline in seqdict[key]:
                file.writelines(f"{seqline}\n")