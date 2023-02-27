#operon finder using .gbk input
import os
import re

def sensecheck(cds):
    if "complement" in cds[1]:
        operonstrand = "negative"
    elif "complement" not in cds[1]:
        operonstrand = "positive"
    return operonstrand

dataset = input("Specify the prokka dataset you will be using('partial', 'comp' or 'both'), or its filepath:")
if dataset == "partial":
    prokkaset = "Data\\Processed_data\\prokka_out_partial_genomes"
elif dataset == "comp":
    prokkaset = "Data\\Processed_data\\prokka_out_complete_genomes"
elif dataset == "both":
    prokkaset = "Data\\Processed_data\\prokka_out_both_genomes"
else:
    prokkaset = dataset

prokkacollection = os.listdir(prokkaset)

for i in range(len(prokkacollection)):
    gbkfile = str(prokkaset + "\\" + prokkacollection[i] + "\\" + prokkacollection[i] + ".gbk")
    with open(gbkfile, "r") as file:
        operongenes = [] #starts nested list for operon constituent genes in format [[LOCUSTAG1, OPERON1],...]
        tmpcoords = [] #nested list that will store location of genes in bases [[gene1start,gene1end],[gene2start,gene2end]...]
        tmp = []
        operoncounter = 0
        linecounter = 0
        iteratenewoperon = True
        for line in file:
            cdsrnacheck = False
            linecounter = linecounter + 1
            splitline = line.split()
            if "CDS" in splitline[0] and " CDS " in line:
                cdsrnacheck = True
            if "tRNA" in splitline[0] and " tRNA " in line:
                cdsrnacheck = True
            if "rRNA" in splitline[0] and " rRNA " in line:
                cdsrnacheck = True
            if "tmRNA" in splitline[0] and " tmRNA " in line:
                cdsrnacheck = True
            if "mRNA" in splitline[0] and " mRNA " in line:
                cdsrnacheck = True
            if "gene" in splitline[0] and " gene " in line:
                cdsrnacheck = True
            if "Location/Qualifiers" in line:
                iteratenewoperon = True
            if cdsrnacheck == True and iteratenewoperon == False:
                newoperonstrand = sensecheck(splitline)
                if originaloperonstrand != newoperonstrand:
                    iteratenewoperon = True
                coords = re.findall("\d+", splitline[1])
                tmpcoords.append(coords)
                genedistance = int(tmpcoords[len(tmpcoords)-1][0]) - int(tmpcoords[len(tmpcoords)-2][1])
                if genedistance > 50:
                    iteratenewoperon = True
                if iteratenewoperon == False:
                    tmp.append("Operon " + str(operoncounter))
            if cdsrnacheck == True and iteratenewoperon == True:
                operoncounter = operoncounter + 1
                iteratenewoperon = False
                originaloperonstrand = sensecheck(splitline)
                coords = re.findall("\d+", splitline[1])
                tmpcoords.append(coords)
                tmp.append("Operon " + str(operoncounter))
            if "locus_tag" in line:
                tagreader = re.findall("locus_tag=\"(.*)\"", line)
                tmp.extend(tagreader)
                operongenes.append(tmp)
                tmp = []
    if dataset == "comp":
        outputtextfile = str("Data\\Processed_data\\operon_predictions_comp\\Operon_prediction_" + prokkacollection[i] + ".txt")
    elif dataset == "partial":
        outputtextfile = str("Data\\Processed_data\\operon_predictions_partial\\Operon_prediction_" + prokkacollection[i] + ".txt")
    elif dataset == "both":
        outputtextfile = str("Data\\Processed_data\\operon_predictions_both\\Operon_prediction_" + prokkacollection[i] + ".txt")
    with open(outputtextfile, "w") as file2: #this may seem clunky vs just making the above outputtextfile = open(...) then doing for elem in opgenes, outputtextfile.writelines(elem), but it kept missing the items at the end of the list if i did it that way
        for element in operongenes:
            file2.writelines(f"{element}\n") 
    progress = (i/(len(prokkacollection)))*100
    print ("The task is " + str(progress) + " percent complete")