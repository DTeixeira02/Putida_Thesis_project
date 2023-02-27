import os
import csv
import re

def blankdel(cell): #function for removing empty cells - this doesnt delete them from the base file, just the rowdata list, errors occurred when blanks present: easier to fix by removing blanks
    if cell == (""):
        return False
    return True  

def scanigrstr(clustlist): #reads the orientation of the IGR for each "GCA_...gene1_IGR_gene2_orient" string from a list made of said string in format [cluster, genome, gene1, gene2, orient]
    generecord = [] #resets the generecord
    temp = [] #resets the temp record
    if clustlist[4] == "CO_R": #each check for clustlist[4] checks the orientation, if its co reverse, co forward etc. further processes decided by orientation
        temp.append(clustlist[0]) #appends the clusterID from [cluster, genome, gene1, gene2, orient]
        temp.append(clustlist[1]) #appends the genome from [cluster, genome, gene1, gene2, orient]
        temp.append(clustlist[2]) #appends the gene1 from [cluster, genome, gene1, gene2, orient]
        generecord.append(temp) #end up with nested list in format [[cluster,genome,gene]] - i know i dont need nested list if i rework code further, its just more familiar for me when writing the output txtfile
    elif clustlist[4] == "CO_F":
        temp.append(clustlist[0])
        temp.append(clustlist[1])
        temp.append(clustlist[3]) #appends the gene2 from [cluster, genome, gene1, gene2, orient]
        generecord.append(temp)
    elif clustlist[4] == "DP":  #different structure as for doublepromoters, the IGR is associated with both neighbouring genes
        temp.append(clustlist[0])
        temp.append(clustlist[1])
        temp.append(clustlist[2])
        generecord.append(temp) #end up with nested list in format [[cluster,genome,gene1]]
        temp = [] #resets temp so can be used again for the remaining neighbour gene
        temp.append(clustlist[0])
        temp.append(clustlist[1])
        temp.append(clustlist[3])
        generecord.append(temp) #end up with nested list in format [[cluster,genome,gene1], [cluster,genome,gene2]]
    elif clustlist[4] == "DT":
        temp.append(clustlist[0])
        temp.append(clustlist[1])
        temp.append("N/A - IGR in/is double terminator")
        generecord.append(temp)
    else:
        generecord.append("Error - No recognised orientation")
    return generecord #returns generecord so that contents can be written to next line(s) of output txtfile

rowdata = []
dataset = input("Specify location of IGR_presence_absence.csv or select 'comp' or 'both':")

if dataset == "both":
    igrcsv = "Data\\Processed_data\\piggy_out_both_genomes\\IGR_presence_absence.csv"
    outputtextfile = open("Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt", "w")
elif dataset == "comp":
    igrcsv = "Data\\Processed_data\\piggy_out_complete_genomes\\IGR_presence_absence.csv"
    outputtextfile = open("Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_complete_genomes.txt", "w")
else:
    igrcsv = dataset
    outputtextfile = open("Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_unspecified_dataset.txt", "w")

print (os.getcwd())
with open(igrcsv, 'r') as file:
    csvreader = csv.reader(file)
    next(csvreader)
    for row in csvreader:
        rowdata.append(row)

for i in range(len(rowdata)):
    currentcluster = list(filter(blankdel, list(rowdata[i])))
    for x in range(len(currentcluster)-4):
        strvar = str(currentcluster[x+4])
        tmp = []
        splittingdoubles = strvar.split()
        if len(splittingdoubles) > 0:
            for b in range(len(splittingdoubles)):
                tmp.append(currentcluster[0])
                splitdelimiters = splittingdoubles[b].split("_+_+_")
                tmp.extend(splitdelimiters)
                generecord = scanigrstr(tmp)
                tmp = []
                for element in generecord: 
                    outputtextfile.writelines(f"{element}\n")
        else:
            tmp.append(currentcluster[0])
            neighbouringgenes = strvar.split("_+_+_")
            tmp.extend(neighbouringgenes)
            generecord = scanigrstr(tmp)
            for element in generecord: 
                outputtextfile.writelines(f"{element}\n")
    progress = (i/len(rowdata))*100
    print ("The task is " + str(progress) + " percent complete")