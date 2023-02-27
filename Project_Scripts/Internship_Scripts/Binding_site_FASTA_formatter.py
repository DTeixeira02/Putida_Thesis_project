#FASTAconverter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import re
import os

locprompt = input("Where is the folder containing raw nucleotide strings?(State 'default' if using for TF text files):")

if locprompt == "default":
    folderloc = Path("Transcription_Factor_Binding_Sites_collecTF") 
    factorfiles = os.listdir(folderloc)

for item in factorfiles:
    recordlist = []
    tfid = re.findall("(.*)_tf_binding_site.txt", item)[0]
    count = 0
    with open(f"{folderloc}\\{item}", "r") as file:
        rawtfstring = [line.strip() for line in file.readlines() if len(line.strip())]
    for x in rawtfstring:
        count = count+1
        recordlist.append(SeqRecord(Seq(x),id=f"{tfid}_{count}"))
    with open(f"TF_binding_site_fasta\\{tfid}_binding_site_fasta.txt", "w") as file:
        for record in recordlist:
            SeqIO.write(record, file, "fasta")