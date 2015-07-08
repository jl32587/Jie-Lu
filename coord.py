#!/usr/bin/python
import re
def create_sra_matrix_from_file_by_len(filename,chr_num,sra_len):
    # accession       start   end     strand  read_length     read_count      hit
    # Chr1    1       23      -       23      1       315

    infile = open(filename)
    infile.readline()
    len = 0

    for line in infile:
        s_line = line.strip()
        b = s_line.split()
        if b[0] == chr_num and b[4] == str(sra_len):
            len += 1
    
    sra_matrix = [[0]*2 for i in range(len)]
    sra_read = [0]*len
    infile.close()
    
    infile = open(filename)
    infile.readline()       #if file has header line
    i = 0
    for line in infile:
        s_line = line.strip()
        a = s_line.split()
        if a[0] == chr_num and a[4] == str(sra_len):
            sra_matrix[i][0] = int(a[1]) # start
            sra_matrix[i][1] = int(a[2]) # end
            #sra_matrix[i][2] = a[3]      # strand
            sra_read[i] = int(a[5])      # reads
            
            i += 1
    infile.close()
    return sra_matrix, sra_read



def create_rep_norm_sra_matrix_from_file(filename,chr_num):
    infile = open(filename)
    infile.readline()
    len = 0

    for line in infile:
        s_line = line.strip()
        b = s_line.split()
        if b[0] == chr_num:
            len += 1
    
    sra_matrix = [[0]*2 for i in range(len)]
    sra_read = [0]*len
    infile.close()
    
    infile = open(filename)
    #infile.readline()
    i = 0
    for line in infile:
        s_line = line.strip()
        a = s_line.split()
        if a[0] == chr_num:
            sra_matrix[i][0] = int(a[1])    # start
            sra_matrix[i][1] = int(a[2])    # end
            sra_read[i] = float(a[5])/float(a[6])  # repeat normalized reads
            i += 1
    infile.close()
    return sra_matrix, sra_read

def create_rep_norm_sra_matrix_from_file_by_len(filename,chr_num,sra_len):
    infile = open(filename)
    infile.readline()
    len = 0

    for line in infile:
        s_line = line.strip()
        b = s_line.split()
        if b[0] == chr_num and b[4] == str(sra_len):
            len += 1
    
    sra_matrix = [[0]*2 for i in range(len)]
    sra_read = [0]*len
    infile.close()
    
    infile = open(filename)
    #infile.readline()
    i = 0
    for line in infile:
        s_line = line.strip()
        a = s_line.split()
        if a[0] == chr_num and a[4] == str(sra_len):
            sra_matrix[i][0] = int(a[1])    # start
            sra_matrix[i][1] = int(a[2])    # end
            sra_read[i] = float(a[5])/float(a[6])  # repeat normalized reads
            i += 1
    infile.close()
    return sra_matrix, sra_read


GFF = "TAIR9_GFF3_genes.gff"
def get_coord_given_agi(agi,feature):

    # gff format:
    # accession       source  feature start   end     score   strand  frame   attributes
    # Chr1    TAIR9   gene    3631    5899    .       +       .       ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
    # feature is a list of interested features

    file_in = open(GFF, 'r')
    file_in.readline()
    coord = {}
    for line in file_in:
        line = line.strip()
        s_line = line.split()
        if s_line[2] in feature:
            attr = s_line[8].split(";")
            tair_agi = re.sub('ID=','',attr[0])
            start = int(s_line[3])
            end = int(s_line[4])       
            coord[tair_agi] = [start,end,s_line[0],s_line[6]] #start,end,chromosome,strand
    file_in.close()
    if agi in coord:
        return coord[agi]
    else:
        print agi+":no coordinate found"


def create_matrix_from_gff(filename,chr):
    'give a matrix of [start,end,strand,agi_identifier] from gff file'
# accession       source  feature start   end     score   strand  frame   attributes
# Chr1    TAIR9   transposable_element    11897   11976   .       +       .       ID=AT1TE00010;Name=AT1TE00010;Alias=ATCOPIA24
    file_in = open(filename)
    file_in.readline()
    gene_matrix = []
    for line in file_in:
        line = line.strip()
        s_line = line.split()
        if s_line[0].upper() == chr.upper():
            start = int(s_line[3])
            end = int(s_line[4])
            attr = s_line[8].split(";")
            tair_agi = re.sub('ID=','',attr[0])
            gene_matrix.append([start,end,s_line[6],tair_agi])
    file_in.close()
    return gene_matrix

def create_matrix_from_kdef_gff(filename,chr):
    'give a matrix of [start,end,strand,agi_identifier] from methylation gff file (Gerhing et al 2009)'
# chrom   source  feature start   end     score   strand  frame   attributes
# chr1    Solexa  WT_embryo_kdef.txt      1201    1301    5.33223847616762e-08    +       .       .

    file_in = open(filename)
    file_in.readline()
    gene_matrix = []
    for line in file_in:
        line = line.strip()
        s_line = line.split()
        if s_line[0] == chr.lower():
            start = int(s_line[3])
            end = int(s_line[4])
            score = float(s_line[5])
            gene_matrix.append([start,end,score])
    file_in.close()
    return gene_matrix

def create_matrix_from_dmr(filename,chr):
    'give a matrix of [start,end,strand,agi_identifier] from dmr file'

#Chr     Context Start   End     EM-EN   p (Fisher's exact test)
#chr1    CHH     11699651        11699800        -0.6400 1.20E-07
#chr2    CHH     17249501        17249550        -0.5000 6.52E-10

    file_in = open(filename)
    file_in.readline() #skip header
    dmr = []
    for line in file_in:
        line = line.strip()
        s_line = line.split()
        if s_line[0].upper() == chr.upper():
            start = int(s_line[2])
            end = int(s_line[3])
            dmr.append([start,end])
            dmr.sort()
    file_in.close()
    return dmr
