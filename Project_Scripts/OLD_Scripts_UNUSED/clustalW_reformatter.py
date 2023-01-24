#clustalW_reformatter

import os
import re

fastqloc = "Thesis_project_data\\Data\\Raw_data\\clustalw.fasta"
outputdir = "Thesis_project_data\\Data\\Raw_data"

with open(fastqloc,"r") as file:
    current_tf = "OmpR"
    fasta_header_to_seq_dict = {}
    fasta_string_list = []
    fasta_line_count = 0
    headerfound = False
    fastafound = False
    firstiteration = True
    for line in file:
        if ">" not in line:
            if len(line.split()) == 1:
                formatted_line = line.strip()
                fasta_string_list.append(formatted_line)
            elif len(line.split()) != 1:
                print(f"Somehow sequence string is discontinuous/contains whitespace {line}")
        if ">" in line:
            if firstiteration == False:
                fasta_string = "".join(fasta_string_list)
                fasta_header_to_seq_dict[header] = fasta_string
                header = line.split()[0]
                fasta_string_list = []
            if firstiteration == True:
                firstiteration = False
                header = line.split()[0]
    with open(str(outputdir+"\\"+current_tf+"_binding_site.fasta"), "w") as file2:
        for header in fasta_header_to_seq_dict.keys():
            file2.write(f"{header}\n{fasta_header_to_seq_dict[header]}\n")
