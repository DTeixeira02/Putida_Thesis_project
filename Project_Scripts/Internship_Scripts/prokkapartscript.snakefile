# This snakemake file is used to generate PROKKA annotation for each of the
# Clostridia genomes

IDS, = glob_wildcards("partialonly/{id}.fna")

rule all:
     input: expand("prokkapartial/{id}/{id}.gff", id=IDS)

# Run PROKKA on an input genome
rule prokka:
     input: "partialonly/{id}.fna"
     output: outdir=directory("prokkapartial/{id}"), outgff="prokkapartial/{id}/{id}.gff",
     shell: "prokka --outdir {output.outdir} --prefix {wildcards.id} --force {input}"
