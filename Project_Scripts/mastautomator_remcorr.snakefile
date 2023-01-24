#mastautomator
genomeids, = glob_wildcards("bothigrs/{id}_IGRs.fna")

rule all:
     input: expand("motifpredbothigrs/{id}/meme.html", id=genomeids), expand("bothigrs/{id}_IGRs.fna", id=genomeids), expand("mastout_remcorr/{id}", id=genomeids)

rule mast:
    input: "motifpredbothigrs/{id}/meme.html","bothigrs/{id}_IGRs.fna"
    output: directory("mastout_remcorr/{id}")
    shell: "mast -mev 0.05 -remcorr -oc {output} {input}"
