# NGS data analysis tools for small RNA-seq in plants

Copyright 2012 Jie Lu. All Rights Reserved.

Small RNAs are 20-30nt long nucleic acids that negatively regulate gene expression/translation. In plants, two major classesof endogenous small RNAs play important roles in post-transcriptional regulation and genome stability maintenance. The 21-nt-long small RNAs are mostly microRNAs which guide target-specific mRNA cleavage. The 24-nt-long small RNAs are generated fromtransposable element-rich regions in the genome and mediate heterochromatin formation and maintenance. 

Next generation sequencing (NGS) has enabled comprehensive
profiling of small RNA populations in plant genomes. This collection of small
RNA-seq analysis tools were developed specifically for the output from CASHX (https://omictools.com/cache-assisted-hash-search-with-xor-logic-tool).

Synopsis:    sRNA.py is an ad hoc Python module for analyzing small RNA-seq NGS data (i.e.18-30nt long RNA molecules that usually negatively regulate gene expression).   
This code reads in .chx file (output from a mapper called CASHX that were specificallydesigned for small RNA-seq, something like .sam file), creates a sRNA class for each library. 
Class functions include generating basic statistics, filtering out unwanted regions, normalizes reads based on both library size and repeats, calls de novo small RNA cluster using a heuristic method,computes normalized reads per cluster, converts to other formats (e.g..wig file which can be visualized in genome browser, etc) and a number of other functions that are customized for different purposes of the study.

Dependency:This code was developed under Python 2.7.2 and is dependent on ol.py, a Python module for finding overlapping intervals using a merge-sort algorithm.

Example: TAIR10_GFF3_genes_transposons.gff is the annotation file of Arabidopsis genome.Usage:A TestRun.py is included to test the codes. This code reads in all .chx files under current directory (a single case "test_input.chx" is provided here), filters out small RNA reads that are mapped to non-small RNA regions (i.e.chloroplast, mitochondria,tRNAs,snoRNAs,snRNAs, and rRNAs), normalizes reads from TEs and repeats, sorts and prints the reads to a .txt file for further analysis.How to Run:python ./TestRun.py

