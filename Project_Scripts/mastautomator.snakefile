#mastautomator
genomeids, = glob_wildcards("bothigrs/{id}_IGRs.fna")

rule all:
     input: expand("motifpredbothigrs/{id}/meme.html", id=genomeids), expand("bothigrs/{id}_IGRs.fna", id=genomeids), expand("mastout/{id}", id=genomeids)

rule mast:
    input: "motifpredbothigrs/{id}/meme.html","bothigrs/{id}_IGRs.fna"
    output: directory("mastout/{id}")
    shell: "mast -oc {output} {input}"
