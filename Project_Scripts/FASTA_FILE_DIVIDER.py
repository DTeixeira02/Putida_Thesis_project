#FASTA_FILE_DIVIDER
import os
import re
from collections import defaultdict
import random

fastqloc = "Thesis_project_data\\Data\\Raw_data\\Unprocessed_SRA_FASTA"
sra_id_to_tf_loc = "Thesis_project_data\\Data\\Raw_data\\SRA_accessions_to_TF.txt"
outputdir = "Thesis_project_data\\Data\\Raw_data\\Processed_SRA_FASTA"
if os.path.isdir(outputdir) == False:
    os.mkdir(outputdir)

with open(sra_id_to_tf_loc,"r") as file:
    id_to_tf_dict = {}
    for line in file:
        splitline = line.split()
        if splitline[0] in id_to_tf_dict.keys():
            print(f"{splitline[0]} appears twice in accession list")
        if splitline[0] not in id_to_tf_dict.keys():
            id_to_tf_dict[splitline[0]] = splitline[1]

ids_done = []
for fasta in os.listdir(fastqloc):
    id = re.findall("(.*)_\d*.fasta",fasta)[0]
    nuc_list = ["A","T","C","G"]
    if id not in ids_done:
        with open(str(fastqloc+"\\"+fasta),"r") as file:
            current_tf = id_to_tf_dict[id]
            fasta_header_to_seq_dict = {}
            fasta_line_count = 0
            headerfound = False
            fastafound = False
            for line in file:
                if ">" not in line:
                    if len(line.split()) == 1 and headerfound == True:
                        if "N" in line:
                            nuc = random.choices(nuc_list,k=1)[0]
                            formatted_line = line.replace("N", nuc)
                            fasta_header_to_seq_dict[header] = formatted_line
                        elif "N" not in line:
                            fasta_header_to_seq_dict[header] = line
                        headerfound = False
                        fasta_line_count = fasta_line_count + 1
                    elif headerfound == False:
                        print(f"Somehow found FASTA sequence before header in {line}")
                    elif len(line.split()) != 1:
                        print(f"Somehow sequence string is discontinuous/contains whitespace {line}")
                if ">" in line:
                    header = line.split()[0]
                    headerfound = True
                if fasta_line_count >= 50:
                    break
            #with open(str(outputdir+"\\"+current_tf+"_binding_site.fasta"), "w") as file2:
            #    for header in fasta_header_to_seq_dict.keys():
            #        file2.write(f"{header}\n{fasta_header_to_seq_dict[header]}")
        ids_done.append(id)
    elif id in ids_done:
        print(f"Already made processed fasta for {id_to_tf_dict[id]}")