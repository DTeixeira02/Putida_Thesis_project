#Locustag_genename_IGR_associator
import os
import re
from collections import defaultdict

def op_pred_reader(operonpred):
    tag_operon_dict = {}
    operon_tag_dict = defaultdict(list)
    with open(operonpred, "r") as file:
        for line in file:
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")
            currentline = currentline.split(",")
            tag_operon_dict[currentline[1]] = currentline[0]
            operon_tag_dict[currentline[0]].append(currentline[1])
    return tag_operon_dict, operon_tag_dict

def igr_pred_reader(igrdoc):
    tag_to_igr_dict = {}
    with open(igrdoc, "r") as file:
        for line in file:
            currentline = line.strip().replace(" ","").replace("'","").replace("[","").replace("]","")  
            currentline = currentline.split(",")
            if len(currentline) == 3 and currentline[2] != "N/A - IGR in/is double terminator":  
                tag_to_igr_dict[currentline[2]] = currentline[0]
            elif len(currentline) != 3:  
                print(f"Potential data issue, {line=}")
    return tag_to_igr_dict

def gbkreader(currentgbk):
    genenamefound = False
    cdsrnacheck = False
    locustagfound = False
    firstiteration = True
    tagtogenedict = {}
    with open(currentgbk, "r") as file:
        for line in file:
            splitline = line.split()
            cdsrnacheck = True
            if "CDS" in splitline[0] and " CDS " in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = False
            if "tRNA" in splitline[0] and " tRNA " in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = False
            if "rRNA" in splitline[0] and " rRNA " in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = False
            if "tmRNA" in splitline[0] and " tmRNA " in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = False
            if "mRNA" in splitline[0] and " mRNA " in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = False
            if "gene" in splitline[0] and " gene " in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = False
            if "ORIGIN" in line:
                if locustagfound == True and genenamefound == True:
                    tagtogenedict[locustag] = genename
                if locustagfound == True and genenamefound == False:
                    tagtogenedict[locustag] = "Unknown_Gene"
                if locustagfound == False and firstiteration == False:
                    print(f"An error has occurred at line {line}, a CDS region lacks a locus tag")
                if genenamefound == True:
                    del genename
                    genenamefound = False
                if locustagfound == True:
                    del locustag
                    locustagfound = False
                firstiteration = True #Intentional, prevents error messaging in fragmented genome datasets
            if "/gene=" in line:
                scangeneline = re.findall("/gene=\"(.*)\"",line)[0]
                checkforwhitespace = scangeneline.split()
                if len(checkforwhitespace) != 1:
                    genename = "_".join(checkforwhitespace)
                elif len(checkforwhitespace) == 1:
                    genename = scangeneline
                genenamefound = True
            if "/locus_tag=" in line:
                locustag = re.findall("/locus_tag=\"(.*)\"",line)[0]
                locustagfound = True
    return tagtogenedict

validchoice = False
while validchoice == False:
    datasetprompt = input("What dataset are you wanting to list genes for?('both' or 'comp'):")
    if datasetprompt == "both" or datasetprompt == "Both" or datasetprompt == "BOTH":
        gbkloc = "Internship_data\\Data\\Processed_data\\prokka_out_both_genomes"
        oploc = "Internship_data\\Data\\Processed_data\\operon_predictions_both"
        igrpredloc = "Internship_data\\Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt"
        outputloc = "Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_both_genomes"
        validchoice = True
    elif datasetprompt == "comp" or datasetprompt == "Comp" or datasetprompt == "COMP":
        gbkloc = "Internship_data\\Data\\Processed_data\\prokka_out_complete_genomes"
        oploc = "Internship_data\\Data\\Processed_data\\operon_predictions_comp"
        igrpredloc = "Internship_data\\Data\\Processed_data\\IGR_gene_associations\\IGR_gene_assocs_both_genomes.txt" #both used intentionally to ensure consistent clusterID between both and comp datasets
        outputloc ="Thesis_project_data\\Data\\Raw_data\\Locustag_to_gene_lists_comp_genomes"
        validchoice = True
    else:
        print("Invalid choice, please choose both/Both/BOTH for 'Both genomes dataset' or comp/Comp/COMP for 'complete genomes dataset'")
genomeids = os.listdir(gbkloc)
tag_to_igr_dict = igr_pred_reader(igrpredloc)

loopcounter = 0
denominator = len(genomeids)
for id in genomeids:
    print(id)
    print(f"Task completion at {(loopcounter/denominator)*100}%")
    valuelist = []
    multiocc = []
    entries_to_write = [] #nested list formatted as [Operon, Locus tag, gene name, clusterID]
    currentgbk = str(gbkloc+"\\"+id+"\\"+id+".gbk")
    currentoppred = str(oploc+"\\Operon_prediction_"+id+".txt")
    taggenedict = gbkreader(currentgbk)
    tag_op_dict, op_tag_dict = op_pred_reader(currentoppred)
    for operon in op_tag_dict.keys(): #Loops for every predicted operon, associates IGR to operon if first operon gene found associated in presenceabsence.csv, else trys last gene, else defines clusterID as NO CLUSTER
        if op_tag_dict[operon][0] in tag_to_igr_dict.keys(): #If first gene in operon is in tag_to_igr_dict keys(dict keyed by gene/locustag)
            clusterID = tag_to_igr_dict[op_tag_dict[operon][0]]
        elif op_tag_dict[operon][len(op_tag_dict[operon])-1] in tag_to_igr_dict.keys(): #If last gene in operon is in tag_to_igr_dict keys(dict keyed by gene/locustag)
            clusterID = tag_to_igr_dict[op_tag_dict[operon][len(op_tag_dict[operon])-1]]
        else:
            clusterID = "NO_CLUSTER_ASSOCIATION_FOUND" #If neither first or last genes in operon associated to IGR
        for locustag in op_tag_dict[operon]:
            entries_to_write.append([operon,locustag,taggenedict[locustag],clusterID])
    totalvalues = list(taggenedict.values())
    for value in totalvalues:
        if value not in valuelist:
            valuelist.append(value)
    for value in valuelist:
        occurrencecount = totalvalues.count(value)
        if occurrencecount > 1:
            print(f"{occurrencecount} occurrences of {value}")
            multiocc.append([value,occurrencecount])
    
    with open(str(outputloc+"\\GENELIST_"+id+".txt"), "w") as file:
        file.write("GENES WITH MULTIPLE OCCURRENCES(Gene -> Occurrence)\n")
        for element in multiocc:
            file.write(f"{element[0]}\t{element[1]}\n")
        file.write("*****************************************************\n")
        file.write("LOCUS TAG TO GENE LIST\n")
        for entry in entries_to_write:
            file.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\t{entry[3]}\n")
    loopcounter = loopcounter+1