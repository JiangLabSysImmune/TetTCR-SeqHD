sample_r1=$1
sample_r2=$2
core=$3

cutadapt -G GTTGTCAAGATGCTACCGTTCAGAG -o tmp/sampletag.R1.fastq.gz -p tmp/sampletag.R2.fastq.gz --pair-filter=both --discard-untrimmed $sample_r1 $sample_r2

if [ -f sampletag_statistics.csv ]; then
        rm -f sampletag_statistics.csv
fi

dir=$(pwd)
cnt1=$(zcat tmp/sampletag.R1.fastq.gz|wc -l)
echo "total sampletag reads,$((cnt1/4))" >> sampletag_statistics.csv
umis fastqtransform --cores $core --separate_cb ref/transform.json tmp/sampletag.R1.fastq.gz tmp/sampletag.R2.fastq.gz|gzip > tmp/sampletag_R2.transform.fastq.gz

umis cb_filter --nedit 2 --cores $core --bc1 ref/bc1.txt --bc2 ref/bc2.txt --bc3 ref/bc3.txt tmp/sampletag_R2.transform.fastq.gz|gzip > tmp/sampletag_R2.filter.fastq.gz

zcat tmp/sampletag_R2.filter.fastq.gz |grep --no-group-separator -A 3 -Ff truecell.txt |gzip > fq/sampletag_R2.true.fastq.gz

cnt2=$(zcat fq/sampletag_R2.true.fastq.gz|wc -l)
echo "total sampletag truecell reads,$((cnt2/4))" >> sampletag_statistics.csv

mkdir sampletag
bowtie2 -p $core --norc --local --no-unal -x $dir/ref/sampletag/sampletag -U fq/sampletag_R2.true.fastq.gz -S sampletag/sampletag.sam

samtools sort -O bam sampletag/sampletag.sam > sampletag/sampletag.sorted.bam
samtools index sampletag/sampletag.sorted.bam
umi_tools count -I sampletag/sampletag.sorted.bam --extract-umi-method umis --per-contig --per-cell --wide-format-cell-counts -S sampletag/sampletag.umi.tsv
