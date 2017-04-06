import sRNA as s
import os
import glob

path = './'
feature = ["protein","transposable_element_gene","transposable_element","transposon_fragment","miRNA","pseudogene"]

for infile in glob.glob( os.path.join(path, '*.chx') ):
    infile = infile.split('/')[1]
    print "Current file is: " + infile
    base = infile.split('.')[0]
    print "Reading cashx file ..."
    s.sRNA.load_cashx(infile)
    print "Filtering chloroplast, mitochondria,tRNAs,snoRNAs,snRNAs, and rRNAs..."
    s.sRNA.remove(mito_chl = "True",feature = ['tRNA','snoRNA','snRNA','rRNA'])
    print "Generating .txt file..."
    s.sRNA.save_txt(base+".txt")
    print "txt file done"
    count = {}
    for i in feature:
        print "Calculating sRNAs that match", i
        total_read = 0
        loci = s.sRNA.get_feature_coord([i])
        count = s.sRNA.count(loci,1,distinct = False)
        for c in count:
            total_read += sum(count[c])
        print total_read

    print "###################################################"


