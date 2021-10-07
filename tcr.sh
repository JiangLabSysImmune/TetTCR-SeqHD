read1=$1
read2=$2
chain=$3
core=$4


if [ -f "$chain"_statistics.csv ]; then
        rm -f "$chain"_statistics.csv
	fi

	cnt1=$(zcat $read1|wc -l)
	echo "total TCR$chain reads,$((cnt1/4))" >> TCR"$chain"_statistics.csv

	umis fastqtransform --cores $core --separate_cb ref/transform.json $read1 $read2|gzip > tmp/TCR"$chain"_R2.transform.fastq.gz

	umis cb_filter --nedit 2 --cores $core --bc1 ref/bc1.txt --bc2 ref/bc2.txt --bc3 ref/bc3.txt tmp/TCR"$chain"_R2.transform.fastq.gz|gzip > tmp/TCR"$chain"_R2.filter.fastq.gz

	zcat tmp/TCR"$chain"_R2.filter.fastq.gz |grep --no-group-separator -A 3 -Ff truecell.txt > fq/TCR"$chain"_R2.true.fastq
	cnt2=$(wc -l < fq/TCR"$chain"_R2.true.fastq)
	echo "total TCR$chain truecell reads,$((cnt2/4))" >> TCR"$chain"_statistics.csv

	cp truecell.txt tcr$chain/
	cd tcr$chain
	bash run.sh ../fq/TCR"$chain"_R2.true.fastq $chain $core > log
	cd ../
	python src/tcr_index_conversion.py tcr$chain/filter.tab out/tcr"$chain"_final.csv

