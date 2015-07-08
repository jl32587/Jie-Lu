#!/usr/bin/python
import re

def print_agi_from_gff(gff_file):
    # gff format:
    # accession       source  feature start   end     score   strand  frame   at
#tributes
    # Chr1    TAIR9   gene    3631    5899    .       +       .       ID=AT1G010
#10;Note=protein_coding_gene;Name=AT1G01010

    file_in = open(gff_file)
    outname_line = gff_file.split(".")
    outname = outname_line[0]+"_agi.txt"
    file_out = open(outname,"w")
    file_in.readline()
    for line in file_in:
        line = line.strip()
        s_line = line.split()
        attr = s_line[8].split(";")
        tair_agi = re.sub('ID=','',attr[0])
        print >>file_out,tair_agi,"\n",

    file_in.close()
    file_out.close()

#print_agi_from_gff("gene_w_te_2k.txt")
#print_agi_from_gff("gene_w_te_body.txt")

def get_strand_given_agi(agi_list,gff_file):
    "get the strand for the agi list and return a dictionary"
    file_in = open(gff_file)
    file_in.readline()
    all_strand = {}
    given_agi_strand = {}
    for line in file_in:
        line = line.strip()
        s_line = line.split()
        attr = s_line[8].split(";")
        tair_agi = re.sub('ID=','',attr[0])
        all_strand[tair_agi] = s_line[6]

    for agi in agi_list:
        given_agi_strand[agi] = all_strand[agi]

    return given_agi_strand

    
def get_gff_from_agi(agi_file,gff_file,outname):
    f1 = open(gff_file)
    f2 = open(agi_file)
    
    out = open(outname,"w")
    f1.readline()
    agi = {}
    for line in f2:
        line = line.strip()
        
        if line:            
            line = line.upper()
            print line
            agi[line] = ""

    for line in f1:
        s_line = line.strip()
        s_line = s_line.split()
        attr = s_line[8].split(";")
        tair_agi = re.sub('ID=','',attr[0])
        if tair_agi in agi:
            print >>out,line, 
