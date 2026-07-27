#!/bin/bash
# bash GRIDseq_analysis.sh WT1.fastq.gz
# Yujie Liu, 6/25/2022

# get parameters, or can set here
fq=$1
name=${1%.fastq.gz}
ncpus=8
bowtie_index=/Share/home/qiyijun/00.project/01.lyj/00.DBs/TAIR10/bowtie_index/TAIR10
chr_len=/Share/home/qiyijun/00.project/01.lyj/00.DBs/TAIR10/TAIR10.chrom_sizes

# remove PCR duplication and adapter
zcat $fq | ./fastx_toolkit_0.0.13/fastx_collapser -Q 33 | ./fastx_toolkit_0.0.13/fastx_clipper -a AGATCGGAAGAGCACACGTCT >$name.clean.fa

# tag and remove linkers, filter read length
python GRIDseq_splitRNAandDNA.py $name.clean.fa

# map to genome, convert to bam (RNA and DNA-end)
for type in RNA DNA; do
    bowtie -p $ncpus -v 0 -m 1 --no-unal -x ${bowtie_index} -f ${name}_${type}.fa -S ${name}_${type}.sam
    samtools view -@ $ncpus -bh -o ${name}_${type}.bam ${name}_${type}.sam
    samtools sort -@ $ncpus -o ${name}_${type}.sort.bam ${name}_${type}.bam
    samtools index -@ $ncpus ${name}_${type}.sort.bam
    rm -f ${name}_${type}.bam
done

# build interactions (RNA to DNA)
bedtools bamtobed -i ${name}_RNA.sort.bam >${name}_RNA.bed
bedtools bamtobed -i ${name}_DNA.sort.bam >${name}_DNA.bed

cat ${name}_RNA.bed ${name}_DNA.bed | awk 'BEGIN{OFS="\t"}; {if(a[$4]){print a[$4],$0}else{a[$4]=$0}}' | awk 'BEGIN{OFS="\t"}; {x=$1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9; a[x]=a[x]?a[x]","$4:$4; b[x]=$6"\t"$12}; END{for(x in a){print x,a[x],".",b[x]}}' >$name.bedpe

echo $name DONE
exit 0
