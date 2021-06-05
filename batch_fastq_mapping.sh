#!/bin/bash

#first processing

module load bwa
module load samtools
module load picard
module load python
module load pylauncher
module load gatk

pfx=$1
splitSize=$2
batchSize=$3
coresPerRun=$4
refDir="/work/00001/mattcowp/Hs_reference_datasets"
#hgReference="$refDir/Homo_sapiens.GRCh37.72.dna.fa"
hgReference="$refDir/human_g1k_v37.fasta"
queue="normal"

echo "-----"r
echo "$pfx:" `date` > $pfx.log


# split and clean up the fastq files                                                                                                                                
rm split.*
echo "" >> ../$pfx.log
echo "python /scratch/02568/mshpak/Share/MERGED_FASTQ/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
python /scratch/02568/mshpak/Share/MERGED_FASTQ/launcher.py split.script 2
echo "Fastq split finished :" `date` >> ../$pfx.log
rm -rf pylauncher_tmp*


#  check to ensure that we go an equal number of splits, if there
# are no r1 or r2 files, then try to re-run split script
numR1Files=`ls r1.* | wc -l`
echo "Number of r1 files: " $numR1Files  `date` >> ../$pfx.log
numR2Files=`ls r2.* | wc -l`
echo "Number of r2 files: " $numR2Files  `date` >> ../$pfx.log
if [ $numR1Files -eq 0 -o $numR2Files -eq 0 ]
then
    echo "python /scratch/02568/mshpak/Share/MERGED_FASTQ/split_fastq_threaded.py $pfx $splitSize"  `date` > split.script
    echo "Initial split seems suspect:  r1=$numR1Files, r2=$numR2Files"  `date` >> $pfx.log
    echo "Waiting while splitting and cleaning the fastq file."  `date` >> ../$pfx.log
    python /scratch/02568/mshpak/Share/MERGED_FASTQ/launcher.py split.script 2
    echo "Fastq split finished :" `date` >> ../$pfx.log
    rm -rf pylauncher_tmp*
    numR1Files=`ls r1.* | wc -l`
    echo "Number of r1 files: " $numR1Files  `date` >> ../$pfx.log
    numR2Files=`ls r2.* | wc -l`
    echo "Number of r2 files: " $numR2Files  `date` >> ../$pfx.log
fi

# check for equal numbers of files
if [ $numR1Files -ne $numR2Files ]
then
    echo "Something wrong with split or fastq files.  Unqueal number of splits."  `date` >> ../$pfx.log
    exit 1;
fi

# 2. Move a sets of splits into sub-directories to be processed independently
echo "Dividing splits into $batchSize per batch: " `date` >> ../$pfx.log
batchNum=0
subdir=""
for i in `seq 0 $((numR1Files-1))`; do
    if [ `expr $i % $batchSize` -eq 0 ]
    then
        subdir="b.$batchNum"
        mkdir $subdir
        batchNum=$((batchNum + 1))
    fi
    mv r2.$i $subdir
    mv r1.$i $subdir
done

# And the residual set, if any:
for file in $( ls r1.* ) ; do
    subdir="b.$batchNum"
    mkdir $subdir
    mv r1.* b.$i
    mv r2.* b.$i
done

# 3. Launch step1_bwa_alignment.sh on each split within its own directory; store job numbers; launch exome_step2.sh to combine chr files
#input arguments for step1: pfx (filename of fastq), $fileExt (extension), $hgreference (reference genome file), $subdir (working directory name)
subdirList=$( ls -d b.* )
rm map.script
for subdir in $subdirList ; do
  echo "Creating launcher for all files in $subdir: " `date` >> ../$pfx.log
  cd $subdir
  for file in $( ls r1.* ); do
     fileExt="${file##*.}"
     echo "Run max_exome_step1.sh on r1.$fileExt and r2.$fileExt" >> ../$pfx.log
     echo "/scratch/02568/mshpak/Share/MERGED_FASTQ/step1_bwa_alignment.sh $pfx $fileExt $hgReference $subdir >& $subdir/mapped.$fileExt.log" >> ../map.script
  done
  cd ..
done

echo "Waiting for set of mapping runs to finish: " `date` >> ../$pfx.log
python /scratch/02568/mshpak/Share/MERGED_FASTQ/launcher.py map.script $coresPerRun
rm -rf pylauncher_tmp*
echo "Mapping Finished: " `date` >> ../$pfx.log

# 4. Launch job to combine final chr files across all directories
# split bam by chromosome
rm -f merge.script
chrList=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
for c in ${chrList[@]};
do
   echo "Merging chr$c "
   echo "samtools merge -f $c.bam b.*/$c.mapped.*.sorted.bam; samtools sort $c.bam $c.sorted; samtools index chr$c.sorted.bam" >> ./merge.script
done
echo "Waiting for per-chromosome merge to finish: " `date` >> ../$pfx.log
python /scratch/02568/mshpak/Share/MERGED_FASTQ/launcher.py merge.script $coresPerRun
rm -rf pylauncher_tmp*
echo "Chromosome merge Finished: " `date` >> ../$pfx.log

# 5. Launch GATK on each chromosome, sorted bam file to realign around indels and do BQSR
rm -f variants.script
for c in ${chrList[@]};
do
   echo "GATK via max_exome_step2.bash on chr$c.sorted.bam"
   echo "/work/02568/mshpak/NEXT_GEN/Testing/max_exome_step2.bash chr$c.sorted chr$c >& variants.chr$c.log" >> variants.script
done
echo "Waiting for variant calling to finish"
python /work/02568/mshpak/NEXT_GEN/Testing/launcher.py variants.script $coresPerRun
rm -rf pylauncher_tmp*
echo "Variant Calling Finished: `date` "

# 6. Merge all bam files into single sample bam
echo "Waiting for final whole-sample merger to finish: " `date` >> ../$pfx.log
samtools merge -f $pfx.sort.bam *.sorted.bam
echo "Whole-sample merge Finished: " `date` >> ../$pfx.log
samtools index $pfx.sort.bam

echo "Processing completed at: " `date` >> ../$pfx.log
exit;
