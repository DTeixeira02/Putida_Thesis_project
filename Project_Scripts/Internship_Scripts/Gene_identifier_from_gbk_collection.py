#Lists all gene names found in all .gbk files
from pathlib import Path
import os
import re

#curline = line.strip().strip("/gene=").strip("\n").replace("\"", "")
def getgenes(currentfile):
    tmpgenelist = []
    with open(currentfile, "r") as file:
        for line in file.readlines():
            if "/gene" in line:
                curline = re.findall("/gene=\"(.*)\"", line)[0]
                tmpgenelist.append(curline)
    return tmpgenelist

locationprompt = input("Where is file containing the prokka output files?:")
if locationprompt == "prokkaboth":
    fileloc = Path("GENOMES/prokkaboth")
prokkafiles = os.listdir(fileloc)
generecord = []
for x in prokkafiles:
    currentgbk = Path.joinpath(fileloc, x, f"{x}.gbk")
    tmpgenes = getgenes(currentgbk)
    for y in tmpgenes:
        if y not in generecord:
            generecord.append(y)
    print(f"Currently on {x}")
print(*generecord)