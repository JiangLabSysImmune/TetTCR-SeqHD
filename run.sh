rhapsody=$1
tetR1=$2
tetR2=$3
tcraR1=$4
tcraR2=$5
tcrbR1=$6
tcrbR2=$7
sampleR1=$8
sampleR2=$9
python src/barcode_conversion.py "$rhapsody"
mkdir fq
mkdir tmp
mkdir out
mkdir figures
cp peptide.csv out/
nohup bash tcr.sh "$tcraR1" "$tcraR2" a 15 &
nohup bash tcr.sh "$tcrbR1" "$tcrbR2" b 15 &
nohup bash sample.sh "$sampleR1" "$sampleR2" 10 &
bash tet.sh "$tetR1" "$tetR2" 10 peptide.csv
umi_tools count -I tet/tetramer_align.dbec.sorted.bam --extract-umi-method umis --gene-tag=GX --per-cell --method unique --wide-format-cell-counts -S tet/tetramer_dbec.tsv
python src/umitools_index_conversion.py tet/tetramer_dbec.tsv out/tetramer_dbec.csv
python src/umitools_index_conversion.py sampletag/sampletag.umi.tsv out/sampletag_umi.csv
source activate scanpy
python src/gene_expression.py "$rhapsody"
source deactivate scanpy
wait
Rscript src/sampletag.R ./out/sampletag_umi.csv "$rhapsody"
Rscript src/tetramer_run.R "$rhapsody"
