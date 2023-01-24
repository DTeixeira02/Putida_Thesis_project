IGR_assoc_analysis contains the outputs of the cluster-gene_association_comparator_V3 and IGR_gene_association_variation_analysis scripts. Information within includes:
 - Gene to IGR association statistics across all genomes
 - Gene to IGR association data PER genome
 - TF-binding-IGR to associated Operon composition data
 - A list of genes that associate with more than one IGR that bind *different* transcription factors

mast_merged_igrs_bothgenomes/mast_merged_igrs_compgenomes contains 1 file per genome as one of the outputs of motif_comparator_V2, each file lists the IGRs that were merged within that genome on the basis of similar motif composition, IGRs on the same line were merged.

mastout contains the mast output files from each genome after running the corresponding meme.html file and FASTA_IGRs_shortenedheaders from Raw_data through MAST

motifpredbothigrs/motifpredcompigrs contains the MEME output for each genome after running the corresponding FASTA_IGRs_shortenedheaders through MEME





V2/V3 in the remaining file headers refer to if the data originally was based on complete graphs produced by either V2 or V3 of the complete_network_creator_operon_to_igr_{VERSION}, respectively.

{version}_Network_graphs_{genomeset}_mastupdated folders contain the output of motif_comparator_V2(Same version of comparator for V2/3 complete graphs) using mastout file corresponding to genome and "complete" network of corresponding genome as inputs, wherein IGRs have been merged according to motif composition and similarity. 

{version}_Graphml_full_regulatory_networks_{genomeset} contains the partially completed graphs produced by complete_network_creator_operon_to_igr_{version} which takes nhmmer data and the incomplete graphs to form edges from operon to IGR on the basis that operon encodes a TF predicted to bind the IGR.