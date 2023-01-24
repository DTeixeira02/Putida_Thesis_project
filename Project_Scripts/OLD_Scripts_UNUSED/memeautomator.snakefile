#MEME motif finder

genomeids, = glob_wildcards("bothigrs/{id}_IGRs.fna")

rule all:
     input: expand("motifpredbothigrs/{id}/{id}_prediction", id=genomeids)

rule memesuite:
    input: "bothigrs/{id}_IGRs.fna"
    output: outdir=directory("motifpredbothigrs/{id}"), outfile="motifpredbothigrs/{id}/{id}_prediction",
    shell: "meme -oc {output.outdir} -dna -mod anr -nmotifs 50 {input}"
