# bla_hmm
Scanning of prokaryotic genomes for beta-lactamase orthologous genes
Scripts for monitoring beta-lactamase spread in prokaryotic genomes using Hidden Markov Models
See the following paper for details:


Operating sytem
CentOS 7

Required software
Python 2.7
HHMER v3.1
CD-HIT v4.8
PROKKA v1.1.2
Blast+ v2.9
SignalP v5

USAGE INSTRUCTIONS

1. BLA_HMM: Identification of beta-lactamase orthologous genes in prokaryotic genomes

Copy the bla_hmm folder to your home directory

the program takes as input fasta files and outputs detected beta-lactamases from the known molecular classes.

Usage:

csh bla_hmm_v11.csh fasta_file beta-lactamase_class(all, A, C, D, B1, B2 or B3; default: all) strain_name(default: strain) output_folder(default: bla_hmm_results) translation_frames(default: 6)

example: csh bla_hmm_v11.csh genome.fasta all

In case genome annotation is available include the *.faa and *.gff files in the working folder. Otherwise the genome will be annotated by PROKKA

2. Genome Scanning:
Scripts for detecting beta-lactamases in NCBIs genomes using bla_hmm

Copy the contents of the genome_scanning folder into your working directory

Retrieve NCBI current taxonomy data for prokaryotic organisms by running the get_fulllineagenames.csh script in the taxonomy folder

Download NCBIs assembly summary files for Archaea or Bacteria using the ftp site

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt --output-document=archaea_assembly_summary.txt

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt --output-document=bacteria_assembly_summary.txt

Run the master_script.csh as follows

csh master_script.csh domain_name(bacteria, archaea or prokaryotes) assembly_summary_file #cpus
Minimum number of cpus required is 8 but use as many as possible

example
csh master_script.csh bacteria bacteria_assembly_summary.txt 8

Results can be found in the following files

"domain"_positive_total.dat
"domain"_greys_total.dat

Only the hits in the positive file can be analyzed with safety.
Hits in the "greys" file require experimental verification.

3. Post Analysis:
Scripts for performing clustering, assessment of catalytic machinery integrity and signal peptide predictions on the identified beta-lactamase hits (all molecular classes) from genome scanning.

Copy the contents of post_analysis folder into your working directory where the results files from genome scanning are located

usage
csh post_analysis.csh results_file(e.g. bacteria_positive_total.dat) #cpus

The script creates one folder per beta-lactamase class named class_"bla_class"_"results_file"
The file containing the data on intact open reading frames is named as: "bla_class"_whole_proteins_positives_meta.dat

4. Hits metrics:
Scripts for measuring the beta-lactamase frequency and richness to the various prokaryotic operational taxonomic units (OTUs) of ranked lineages.

Copy the contents of the hits_metrics folder into your working directory containing the results file from genome scanning and the post analysis folders for each beta-lactamase class

usage
csh hits_statistics_v2_ranked.csh "results_file" "output_folder_name"

Beware that the script creates additional results files where the unranked OTUs have been removed from the respective lineages.

5. Protein cluster metrics:
Scripts for assessing the organisms to which each protein cluster belongs to as well as their secretion potential and catalytic machinery integrity.

Copy the contents of the protein_cluster_metrics folder into your working directory containing the data from the previous analyses.

First run the protein_clusters_taxa_multi.csh as follows
csh protein_clusters_taxa_multi.csh "bla_class" "results_file" #cpus
e.g. csh protein_clusters_taxa_multi.csh A bacteria_positive_total.dat 8
Then run the proteins_clusters_data.csh as follows
csh proteins_clusters_data.csh "bla_class" "results_file"

The final script outputs a result file containing the data foreach cluster named data_protein_clusters.dat located in the folder of each beta-lactamase class

6. iTol tree annotation:
Scripts for annotating trees of the inferred phylogenies of the identified protein clusters. Only class A, B1 and B2 scripts are given. In case of other beta-lactamase class the respective lines for assessing the presence of catalytic amino acids should be modified. Files containing HEX color codes for bacterial classes and phyla are also given.

Copy each script along with its color codes file into the respective beta-lactamase class folder containing the results from the previous analyses.

run the script as follows:
csh B1_itol_annotation_exa.csh data_protein_clusters.dat "preffix of each protein cluster used during alignment"

The scripts output three files containing the data that should be added in the respective templates of iTol
1. for_classes_color_strip.dat containing the bacterial class (in case > 1 the clade is annotated as MULTIPLE_CLASSES) of each cluster and the respective color for the color strip template
2. for_multibar.dat containing the size of each cluster (log(#members)+1) and the number of different genera (log(#genera)+1) that should be added in the multibar template
3. for_binary.dat containing binary data on the presence of all catalytic residues, secretion potential and experimental verification of beta-lactam hydrolysis

