import re
import os
import ol
from random import sample

class seq(type):
    def __iter__(self):
        for i in self.sRNA_map:
            yield i

class sRNA(object):
    # sRNA class do not have __init__() method because sRNA dataset is huge and it is not memory efficient to create instances for each sRNA seq
    __metaclass__ = seq
    sRNA_map = []
    norm_base = 1000000
    total_size = 0
    unique_size = 0
    @staticmethod
    # a yield function
    def generator(x):
        for i in x:
            yield i

    @staticmethod
    def get_item(line):
        line = line.strip()
        return line.split()

    @staticmethod
    def sort_map(seq_map):
        '''sort a mapped sequence by its acc, start, end'''
        seq_map.sort(key = lambda x: (x[0], x[1], x[2]))

    @classmethod
    def __iter__(cls):
        for i in cls.sRNA_map:
            yield i
    
    @classmethod
    # load data from cashx output and add it to the sRNA list
    def load_cashx(cls, cashx_file,TAG = False):
        '''load .cashx file into class sRNA and create a list of sRNA and returns library size (total, unique)'''
        #ID=ILLUMINA-755F90_0001:2:48:1280:1888#0/1|reads=1      HITS=52394 SEQ=AAAAAAAAAAAAAAAAAA
        #Accession=NODE_53734_length_1731_cov_3.495090   Start=450       End=468 Strand=1
        #Accession=NODE_53734_length_1731_cov_3.495090   Start=449       End=467 Strand=1
        #Accession=NODE_53734_length_1731_cov_3.495090   Start=448       End=466 Strand=1
        #Accession=NODE_41136_length_698_cov_0.984241    Start=702       End=720 Strand=1
        cls.total_size = 0
        cls.unique_size = 0
        if len(cls.sRNA_map) > 0: del cls.sRNA_map[:]
        f = open(cashx_file, 'r')
        for line in cls.generator(f):
            if line.startswith('ID='): ## line of read information
                item = cls.get_item(line)
                #item[0] = re.sub('ID=ILLUMINA-','', item[0])
                #id = item[0].split('|')[0]
                read = int(re.sub('reads=', '', item[0].split('|')[1]))
                hit = int(re.sub('HITS=', '', item[1]))
################library size is calculated when the library is loaded##########
                cls.total_size += read    ## calculate total lib siz
                cls.unique_size += 1      ## calculate unique lib size
                if TAG == True:
                    try: # if cashx file has SEQ printed out
                        tag = re.sub('SEQ=', '', item[2])
                    except:
                        print "TAG is missing"
                        pass                
            else:         ## line of map information
                item = cls.get_item(line)
                acc = re.sub('Accession=', '', item[0])
                if acc.startswith('NODE'):  ## if database is Aa contigs, only get the node number
                    acc = acc.split('_')[1]
                start = int(re.sub('Start=', '', item[1]))
                end = int(re.sub('End=', '', item[2]))
                strand = int(re.sub('Strand=', '', item[3]))
                if strand == 1:
                    strand = '+'
                else:
                    strand = '-'
                if TAG == True:
                    cls.sRNA_map.append((acc, start, end, read, hit, strand, tag))
                else:
                    cls.sRNA_map.append((acc, start, end, read, hit, strand))
                # sRNA_map index:
                # 0: acc
                # 1: start
                # 2: end
                # 3: read
                # 4: hit
                # 5: strand
                # 6: tag    ## optional
        cls.sort_map(cls.sRNA_map) # sort by acc, start, end
        f.close()

    @classmethod
    def save_gff(cls, filename):
        '''save sRNA list to .gff file'''
        # specification of .gff format for sRNA class:
        # 1. seqname (chromosome or contig number)
        # 2. source -- 'jl'
        # 3. sRNA size
        # 4. start
        # 5. end
        # 6. reads
        # 7. strand
        # 8. frame "."
        # 9. hits
        f_out = open(filename, 'w')
        for i in cls.generator(cls):       
            print >>f_out, '%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d' % (i[0], 'jl', i[2] - i[1] + 1,
                                                                   i[1], i[2], i[3], i[5],
                                                                   '.', i[4])            
        f_out.close()        

    @classmethod
    def save_tab(cls, filename):
        '''save sRNA list to .txt file which can be the tab format'''
        # tab format column:
        # 1.chr
        # 2.start
        # 3.end
        # 4.strand
        # 5.length
        # 6.count
        # 7.hit
        f_out = open(filename, 'w')
        print >>f_out,  "accession\tstart\tend\tstrand\tread_length\tread_count\thit"
        for i in cls. generator(cls):
            print >>f_out, 'Chr%s\t%d\t%d\t%s\t%d\t%d\t%d' % (i[0],i[1],i[2],i[5],i[2]-i[1]+1,i[3],i[4])
        f_out.close()
        
    @classmethod
    def save_segmentSeq(cls, filename):
        '''save sRNA list to .txt file which can be the input of segmentSeq package'''
        # segmentSeq format column:
        # 1.chr
        # 2.seq
        # 3.count
        # 4.start
        # 5.end
        # 6.strand
        f_out = open(filename, 'w')
        print >>f_out,  'chr\ttag\tcount\tstart\tend\tstrand'
        for i in cls. generator(cls):
            print >>f_out, '>Chr%s\t%s\t%d\t%d\t%d\t%s' % (i[0],i[6],i[3],i[1],i[2],i[5])
        f_out.close()

    @classmethod
    
    def save_wig(cls, lib, ref_genome = "At", factor = None, step = 10, unique = True):
        '''to save sRNA in wig format for gbrowse'''
        filename = lib+'_'+ref_genome+'.wig'  # create output filename
        f_out = open(filename, 'w')
        track_line = '''track type=wiggle_0 name="%s" description="%s" visibility=full color=255,0,0 altColor=0,100,200 priority=20''' % (lib, lib)       ## create track line
        if factor == None:# if no factor is given, calculate it from library size
            if unique == True:  ## calculate normalization factor
                factor = float(cls.unique_size)/cls.norm_base
            else:
                factor = float(cls.total_size)/cls.norm_base
        print factor
            
        sRNA_chr = list(set([i[0] for i in cls]))
        sRNA_chr.sort()
        print >>f_out, track_line
        for c in sRNA_chr:
            if ref_genome == "At": # if reference genome is A.thaliana, need to add 'Chr'
                def_line = 'fixedStep chrom=Chr%s start=1 step=%s span=%s' % (c, step, step)
            else:
                def_line = 'fixedStep chrom=%s start=1 step=%s span=%s' % (c, step, step)
            a = filter(lambda x: x[0] == c, cls)
            total_bp = a[-1][2] + 1
            value = [0 for i in range(total_bp+1)]
            bin = total_bp / step + 1
            for s in cls.generator(a):  # for sRNAs on c chromosome
                for i in range(s[1], s[2]+1):
                    if unique == True:
                        try:
                            value[i] += float(1)/s[4]     #unique reads
                        except IndexError:
                            print 'i=%s' % i
                            print s
                    else:
                        value[i] += float(s[3])/s[4]  #total reads
            value = [sum(value[i*step:(i+1)*step]) for i in range(0,bin)]
            value.append(sum(value[(bin-1)*step:]))
            print >>f_out, def_line    
            for n in value:
                print >>f_out, '%0.2f' % (float(n)/factor)
        f_out.close()
        
    @classmethod
    def load_gff(cls, gff):
        cls.sRNA_map = []
        f = open(gff, 'r')
        for line in cls.generator(f):
            item = cls.get_item(line)
            cls.sRNA_map.append((item[0], int(item[3]), int(item[4]),
                                 int(item[5]), int(item[8]), item[6]))
            ## start, end, read, hit should be integers
            cls.sort_map(cls.sRNA_map)
        f.close()

    @classmethod
    def load_tab(cls, tab):
        cls.sRNA_map = []
        f = open(tab, 'r')
        f.readline() # skip header line
        for line in cls.generator(f):
            item = cls.get_item(line)
            cls.sRNA_map.append((item[0], int(item[1]), int(item[2]),
                                 int(item[5]), int(item[6]), item[3]))
        f.close()
                                                          
    @classmethod
    def get_feature_coord(cls,*feature):
        annotation_file = "TAIR9_GFF3_genes.gff"
        f = open(annotation_file, 'r')
        f.readline()  # get rid of the header
        coord = []
        for line in cls.generator(f):
            item = cls.get_item(line)
            if item[2] in feature:
                acc = re.sub('Chr','',item[0])
                start = int(item[3])
                end = int(item[4])
                coord.append((acc, start, end))
        f.close()
        return coord[:]
    
    ##get a list of coordinates given a list of features chosen from following:
    ##chromosome, gene, mRNA, protein, exon, five_prime_UTR, CDS, three_prime_UTR
    ##miRNA, tRNA, ncRNA, snoRNA, snRNA, rRNA
    ##pseudogene, pseudogenic_transcript, pseudogenic_exon
    ##transposable_element_gene,  transposable_element, transposon_fragment
    
    ##gff file specifications for annotation file:
    ##1. accession
    ##2. source
    ##3. feature
    ##4. start
    ##5. end
    ##6. score
    ##7. strand
    ##8. frame
    ##9. attributes

    @classmethod
    def stat(cls):
        '''provide statistics of the sRNA library'''
        total_r = sum([float(i[3])/i[4] for i in cls.generator(cls)])
        unique_r = sum([float(1)/i[4] for i in cls.generator(cls)])
        size_dist = {}
        for l in range(20,26):
            size_dist[l] = [0,0]
        for i in cls:
            size = i[2] - i[1] + 1
            if size in size_dist:
                size_dist[size][0] += float(1)/i[4]   #unique reads
                size_dist[size][1] += float(i[3])/i[4]#total reads                                                                              
        print 'total mapped reads is: %0.1f' % total_r
        print 'unique mapped reads is: %0.1f' % unique_r
        print 'size distribution of unique or total reads:'
        for i in size_dist.keys():
            print '%s\t%0.2f\t%0.2f' % (i, size_dist[i][0]/unique_r, size_dist[i][1]/total_r)
            

    @classmethod
    def find(cls, locus, remove = False):
        '''find all sRNAs mapped to a given genomic region'''
        # locus needs to be sorted before doing overlap
        cls.sort_map(locus)
        locus_chr = list(set([i[0] for i in locus]))
        sRNA_chr = list(set([i[0] for i in cls]))
        locus_chr.sort()
        sRNA_chr.sort()
        new_list = []
        for c in sRNA_chr:
            a = filter(lambda x: x[0] == c, cls)
            del_index = []
            if c not in locus_chr: # if no need to search that chromosome:
                if remove == True: # no sRNA to be removed from sRNA_map   
                    new_list = new_list + a[:]
                next
            elif c in locus_chr:
                b = filter(lambda x: x[0] == c, locus)
                # extract start and end, overlap them
                o = ol.overlap([i[1:3] for i in a],[i[1:3] for i in b])
                if all([len(i) == 0 for i in o]): # if no sRNAs were found to match any given loci
                    new_list = new_list + a[:]
                    print 'no sRNAs were found on Chr%s' % c
                else:
                    for i in xrange(len(o)):  # for each locus
                        sRNA = []
                        try:    # if agi is provided
                            agi = b[i][3]
                            if remove == False:
                                print '@' + agi
                        except IndexError:
                            pass                
                        if o[i]:  # if there is at least one sRNA in that locus
                            del_index = del_index + o[i]
                            for j in o[i]:
                                sRNA.append(a[j])
                            if remove == False:
                                for s in sRNA:
                                    print 'chr%s\t%d\t%d\t%d\t%d\t%s' % s
                        else:     # if no sRNA is mapped to that locus
                            if remove == False:
                                print 'NA'
                    if remove == True: # if sRNAs need to be removed from sRNA_map
                        for k in xrange(len(del_index)-1, -1, -1): # del from the largest index
                            try:
                                del a[del_index[k]]
                            except IndexError:
                                print c, len(a), k
                        new_list = new_list + a[:]
                
        if remove == False:
            for c in locus_chr:
                if c not in sRNA_chr:   ## if acc of a locus is not in sRNA list
                    b = filter(lambda x: x[0] == c, locus)
                    for i in b:
                        try:
                            agi = i[3]
                            print '@' + agi
                            print 'No accession found'
                        except IndexError:
                            pass
            
        elif remove == True:
            cls.sort_map(new_list)
            return new_list[:]
                             
    @classmethod
    def remove(cls, size = None, mito_chl = False, *feature):
        '''to filter the sRNA list with size, chromosome or features,
        note: size is a list of nt that needs to be kept'''
        if size != None:    # size exclusion
            cls.sRNA_map = filter(lambda x: x[2]-x[1]+1 in size, cls.generator(cls))
        if mito_chl == True:   # if mitochondrial and chloroplast sRNAs need to be removed
            cls.sRNA_map = filter(lambda x: x[0].isdigit(),cls.generator(cls))
        if feature:   # if sRNAs from some loci need to be removed
            locus = cls.get_feature_coord(*feature)
            cls.sRNA_map = cls.find(locus, remove = True)

    @classmethod
    def locus_avg_sRNA(cls, locus, flank, window, norm_factor):
        '''calculate the average distinct sRNA density in given loci'''
        cls.sort_map(locus)
        locus_chr = list(set([i[0] for i in locus]))
        sRNA_chr = list(set([i[0] for i in cls]))
        chr = locus_chr & sRNA_chr
        chr.sort()
        locus_num = len(locus)
        point = flank/window
        up = [0]*point
        trans = [0]*point
        down = [0]*point
        
        for c in chr:
            a = filter(lambda x: x[0] == c, cls)
            b = filter(lambda x: x[0] == c, locus)
            o = ol.overlap([i[1:4] for i in a],
                           [[i[1]-flank, i[2]+flank, i[3],i[4]] for i in b])
            for i in range(len(o)):
                if o[i]:
                    locus_len = b[i][2]-b[i][1]+1+flank*2
                    each_locus = [0] * locus_len
                    for j in o[i]:
                        rel_coord = [a[j][1]-b[i][1]+flank,
                                     a[j][2]-b[i][1]+flank]
                        for k in range(max(0,rel_coord[0]),
                                       min(locus_len,rel_coord[1])+1):
                            each_locus[k] += 1
                    if b[i][3] is "-":
                        each_locus.reverse()
                        
                    up_list = each_locus[0:flank]
                    trans_list = each_locus[flank:locus_len-flank]
                    down_list = each_locus[locus_len-flank:]

                    for i in range(point):                   
                        up[i] += float(sum(up_list[i*window:(i+1)*window]))/window
                        down[i] += float(sum(down_list[i*window:(i+1)*window]))/window
                            
                        a = len(trans_list)/point
                        if a >= 1:
                            for i in range(point):
                                trans[i]+= float(sum(trans_list[i*a:(i+1)*a]))/a
        return [i/(locus_num * norm_factor) for i in up + trans + down]

    @classmethod
    def count(cls, loci, norm_factor, distinct = True):
        '''count sRNA reads within given collection of loci
        loci: [chr, start, end] MUST BE SORTED'''
        locus_chr = set([i[0] for i in loci])
        sRNA_chr = set([i[0] for i in cls])
        chr = list(locus_chr & sRNA_chr)
        chr.sort()
        count = {}
        for c in chr:
            count[c] = []
            a = filter(lambda x: x[0] == c, cls)
            b = filter(lambda x: x[0] == c, loci)
            aa = [[i[1],i[2]] for i in a]
            bb = [[i[1],i[2]] for i in b]
            o = ol.overlap(aa, bb)
            for i in range(len(o)):
                if o[i]:
                    read = 0
                    for j in o[i]:
                        if distinct == True:
                            read += float(1)/norm_factor
                        else:
                            read += float(a[j][3])/(a[j][4]*norm_factor)
                    count[c].append(read)
                else:
                    count[c].append(0)
        return count

def cluster(a, max_gap, min_reads, norm_factor):
    '''assign siRNA loci where neighboring sRNA reads are within max_gap and with at least min_reads number of reads and distinct sRNA density within a loci is equal or greater than min_rpkm....a = [[start1, end1], [start2, end2]...[starti,endi]] sorted list'''
    n = 0         # set the initial index of the first sRNA of the first candidate loci
    j = 0         # set the initital index of the current sRNA within the current loci
    reads = 1
    loci_len = 0
    clust_sRNA = []
    while n+j+1 < len(a):
        this = n + j      # keep track of the current sRNA
        next = n + j + 1  # keep track of the next sRNA
        if next < len(a)-1: # if next sRNA is not the end of the list
            if a[next][0] <= a[this][1] + max_gap: # if next sRNA is within max_gap distance
                reads += 1
                j += 1        # go on to the next sRNA
            elif a[next][0] > a[this][1] + max_gap:# if there is no sRNA wthin max_gap right to the current sRNA
                loci_len = a[this][1] - a[n][0] + 1  # calculate sRNA locus length
                if float(reads)/norm_factor >= min_reads: # if the number of NORMALIZED distinct reads within a locus is greater than cutoff
                    clust_sRNA.append([a[n][0], a[this][1]])
                    
                n = next         # start a new candidate loci
                j = 0            # reset the current sRNA within a new candidate locus
                reads = 1        # reset the number of reads within a new locus
                loci_len = 0     # reset the locus length of a new locus
        else: # if next sRNA reaches to the end of the list
            if a[next][0] <= a[this][1] + max_gap:
                reads += 1
                if float(reads)/norm_factor >= min_reads:
                    clust_sRNA.append([a[n][0], a[next][1]])
                    break
                elif float(reads)/norm_factor < min_reads:
                    break    
            elif a[next][0] > a[this][1] + max_gap:
                if float(reads)/norm_factor >= min_reads:
                    clust_sRNA.append([a[n][0], a[this][1]])
                    break
                elif float(reads)/norm_factor < min_reads:
                    break
    return clust_sRNA

def merge_sra_loci(x):
    n = max([i[-1][1] for i in x])
    l = [0] * (n+2)
    m = 0
    merged = []
    for i in x:
        for j in i:
            for k in range(j[0],j[1]+1):
                l[k] += 1
    coord = [0,0]
    if l[0] == 1:   # if the first locus includes the first nucleotide
        coord[0] = 0    
        for i in range(1,n+1):
            if l[i] != 0: 
                if l[i+1] == 0:
                    coord[1] = i
                    merged.append(coord)
                    coord = [0,0]
                else:
                    pass
            elif l[i] == 0:
                if l[i+1] != 0:
                    coord[0] = i+1
                else:
                    pass
    else:         # if the first locus does not include the first nucleotide
        for i in range(1,n+2):
            if l[i] != 0:
                if l[i-1] == 0:
                    coord[0] = i
                else:
                    pass
            elif l[i] == 0:
                if l[i-1] != 0:
                    coord[1] = i-1
                    merged.append(coord)
                    coord = [0,0]
                else:
                    pass
                
    return merged
###############################################################################
## toy dataset for cluster and merge_sra_lci
###############################################################################
## a = [[1,23],[3,25],[4,26],[100,123],[102,125],[226,250],[230,253],[240,263],[450,473],[452,470],[480,503],[1000,1024]]
## print cluster(a,100,3,1)

## a = [[1,23],[100,123],[226,250],[230,253]]

## b = [[2,25],[90,130],[102,125],[240,263],[450,473]]

## c= [[3,30],[89,150],[160,200],[449,800]]

## d = [[0,50]]
## print merge_sra_loci([a,b,c,d])
## print cluster(a,100,3,1)

###############################################################################
## bootstrap to test the enrichment in a certain type of loci
##############################################################################
def z_score(data, x):
    std = lambda d:(sum((x-1.*sum(d)/len(d))**2 for x in d)/(1.*(len(d)-1)))**.5
    avg = lambda x: float(sum(x))/len(x)
    z = (data - avg(x))/std(x)
    return z
    
def resample(n, sra_loci):
    '''n is the length of the chromosome
    sra_loci = [[s1,e1],[s2,e2]...[sn,en]]'''
    loci = []
    start = sample(range(1,n),len(sra_loci))
    start.sort()
    
    for i in range(len(sra_loci)):
        l = sra_loci[i][1] - sra_loci[i][0]
        loci.append([start[i], start[i]+l])
    return loci

def bootstrap(n,chr_len, sra_loci, feature_loci):        
    data = ol.overlap(sra_loci, feature_loci, percentage = True)
    boot = []
    for i in range(n):
        resampled_loci = resample(chr_len, sra_loci)
        ol_perc = ol.overlap(resampled_loci, feature_loci,
                                  percentage = True)
        boot.append(ol_perc)
    z = z_score(data, boot)
    return z

###############################################################################
# toy dataset for bootstrap()
###############################################################################
## chr_len = 500
## n = 100
## sra_loci = [[25,49],[80,99],[150,173],[300,333]]
## genome_loci = [[1,24],[50,73],[100,123],[400,500]]
## print 'z-score is:'
## print bootstrap(n,chr_len,sra_loci,genome_loci)

##############################################################################
# calculate overlap z-score for small RNA loci in mutant crosses
##############################################################################

#sra_f = open('nrpd1a_cross_merged_20.loci','r')
sra_files = ['none_maternal_dependent.loci','none_p4_dependent.loci',
             'strong_maternal_dependent.loci','weak_maternal_dependent.loci']
#loci_f = open('tair9_TE.txt','r')  ## TE
#loci_f = open('gene_w_te_2k.txt','r')    ## TAG
#loci_f = open('protein_coding_gene_minus_te.txt','r') ## non-TAG
#loci_f = open('EP_ESS.gff','r')
#loci_f = open('OST_ESS.gff','r')
#loci_f = open('DMR_CG_EN.txt','r')
#loci_f = open('DMR_CHH_EN.txt','r')
#loci_f = open('DMR_CG_EB.txt','r')
loci_f = open('DMR_CHH_EB.txt','r')

chr_len = [30427671,19698289,23459830,18585056,26975502]
n = 10
loci = []
for line in loci_f:
    line.strip()
    line = line.split()
   # if line[0] == 'Chr1':
    if line[0] == 'chr1':
        
        #loci.append([int(line[3]),int(line[4])])
        loci.append([int(line[2]),int(line[3])])
loci_f.close()

for f in sra_files:
    sra_f = open(f,'r')
    sra_loci = []
    for line in sra_f:
        line.strip()
        line = line.split()
        if line[0] == 'Chr1':
            sra_loci.append([int(line[1]),int(line[2])])

    print 'z-score for %s' % f
    print bootstrap(n,chr_len[0],sra_loci,loci)
    sra_f.close()
    
###############################################################################
#Assign small RNA loci to individual library
##############################################################################

## lib = ['4x4_tag.chx','4xn_tag.chx','nx4_tag.chx','nxn_tag.chx','4xr_tag.chx',##        'rx4_tag.chx','rxr_tag.chx']
## norm_factor = [1.06, 2.16, 1.89, 1.90, 2.25, 1.94, 1,95]

## for i in range(len(lib)):
##     print 'loading...'+lib[i]
##     tab = lib[i].split('.')[0]+'.tab'
##     sRNA.load_cashx(lib[i],TAG=True)
##     print lib[i]+' load done...'
##     sRNA.remove(size = [20,21,22,23,24,25], mito_chl = True)
##     sRNA.save_tab(tab)
##     print "save done"

## for i in range(len(lib)):
##     tab = lib[i].split('.')[0]+'.tab'
##     loci = lib[i].split('.')[0]+'.loci'
##     print "clustering small RNA loci for "+lib[i]
##     norm = norm_factor[i]
##     out = open(loci,'w')
##     sRNA.load_tab(tab)
##     chr = list(set([i[0] for i in sRNA.sRNA_map]))
##     chr.sort()

##     for c in chr:
##         a = filter(lambda x:x[0] == c, sRNA.sRNA_map)
##         a = [[i[1],i[2]] for i in a]
##         b = cluster(a,200,10,norm)
##         print c+' done'
##         for locus in b:
##             print >>out, c, locus[0], locus[1]
##     out.close()


###############################################################################
# merge loci
###############################################################################
## locus_file = ['4x4_tag.loci','4xn_tag.loci','nx4_tag.loci','nxn_tag.loci']
## f_out = open('nrpd1a_cross_merged.loci','w')

## locus_file = ['4x4_tag.loci','4xr_tag.loci','rx4_tag.loci','rxr_tag.loci']
## f_out = open('rdr2_cross_merged_10.loci','w')

## f_in = open(locus_file[0],'r')
## m = {}
## for line in f_in:
##     line.strip()
##     item = line.split()
##     if item[0] not in m:
##         m[item[0]] = []
## print m.keys()

## for f in locus_file:
##     file_in = open(f, 'r')
##     a = []
##     for line in file_in:
##         line.strip()
##         item = line.split()
##         a.append([item[0],int(item[1]),int(item[2])])
##     for c in m.keys():
##         b = filter(lambda x: x[0] == c,a)
##         m[c].append([[i[1],i[2]] for i in b])
## f_in.close()

## chr = m.keys()
## chr.sort()
## for c in chr:
##     merged = merge_sra_loci(m[c])
##     for i in merged:
##         print >>f_out, c, i[0], i[1]
        
## f_out.close()

######################################################################
## calculate the read counts for each sra loci
#####################################################################
# load sra loci
## loci_f = open('nrpd1a_cross_merged_20.loci','r')
## loci = []
## for line in loci_f:
##     item = sRNA.get_item(line)
##     loci.append([item[0], int(item[1]), int(item[2])])
## sRNA.sort_map(loci)
## print 'loci loaded'
## loci_f.close()

## # loci sra data
## lib = {'4x4':1.06,'4xn':2.16,'nx4':1.89,'nxn':1.90}

## for l in lib:
##     ct = {}
##     sRNA.load_tab(l + '_tag.tab')
##     sRNA.remove(size = [24])
##     print 'sra_loaded:' + l
##     print 'norm factor is:' + str(lib[l])
##     out = open(l + '.ct','w')
##     # calculate sra reads for each locus
##     ct = sRNA.count(loci,lib[l])
##     chr = ct.keys()
##     chr.sort()
##     for c in chr:
##         for i in ct[c]:
##             print >>out, c, i
##     out.close()


###############################################################################
# something that might be useful
###############################################################################
## def RPKM(read, norm_factor, loci_len):
##     'calculate RPKM: reads per kilobasepair per TEN million reads' 
##     return float(read)/((float(loci_len)/1000)*norm_factor)

## def get_locus_coord(ref_gff, agi = None, feature = None):
##     f = open(ref_gff, 'r')
##     for line in sRNA.generator(f):
##         item = sRNA.get_item(line)
##         if agi != None:
            
##         locus_coord.append((item[0], int(item[3]), int(item[4]),
##                                  int(item[5]))
##             ## start, end, read, hit should be integers
         
##         f.close()
##     return locus_coord



