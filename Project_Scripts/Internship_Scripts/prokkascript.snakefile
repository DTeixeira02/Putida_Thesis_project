# This snakemake file is used to generate PROKKA annotation for each of the
# Clostridia genomes

IDS, = glob_wildcards("complete/{id}.fna")

rule all:
     input: expand("prokkacomp/{id}/{id}.gff", id=IDS)

# Run PROKKA on an input genome
rule prokka:
     input: "complete/{id}.fna"
     output: outdir=directory("prokkacomp/{id}"), outgff="prokkacomp/{id}/{id}.gff",
     shell: "prokka --outdir {output.outdir} --prefix {wildcards.id} --force {input}"
