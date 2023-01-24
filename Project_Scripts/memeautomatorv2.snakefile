#MEME motif finder

genomeids, = glob_wildcards("bothigrs/{id}_IGRs.fna")

rule all:
     input: expand("motifpredbothigrs/{id}", id=genomeids)

rule memesuite:
    input: "bothigrs/{id}_IGRs.fna"
    output: directory("motifpredbothigrs/{id}")
    shell: "meme -oc {output} -dna -mod anr -nmotifs 50 {input}"
