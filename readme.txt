Required input files:
1.	Fastq files:
Four fastq files should be obtained after sequencing:
i)	Gene expression + AbSeq + SampleTag
ii)	Tetramer barcodes
iii)	TCRa
iv)	TCRb
2.	Peptide reference file:
The peptide reference file should be in csv format and named as peptide.csv. Consisting of five columns (See below). Make sure the column header is the exactly same as below. Different categories are typically labeled on different fluorophores. If there is only one category, it’s OK. Also make sure during experiments, there is Empty barcode(s) as UV-exchange negative control.
You CAN have additional columns after the category column for your custom analysis.
name	Ntsequence	AAsequence	gene	category
IA2	TTCATTGACTCTTACATCTGCCAGGTT	SLSPLQAEL	IA2	PE
HCV	AAACTGGTTGCTCTGGGTATCAACGCTGTT	KLVALGINAV	HCV	APC
Empty	AACCTGGTTCCGATGGTTGCTACCGTT	Empty	Empty	Empty

Required tool installation:
1.	Sequencing processing:
i.	Cutadapt (https://cutadapt.readthedocs.io/en/stable/guide.html), v1.8 (for python 2) or v2.1 (for python 3)
ii.	umis (https://github.com/vals/umis)
iii.	UMI-tools (https://github.com/CGATOxford/UMI-tools)
iv.	Seqtk (https://github.com/lh3/seqtk)
v.	Samtools (http://samtools.sourceforge.net/)
vi.	Scanpy (https://icb-scanpy.readthedocs-hosted.com/en/latest/index.html; NOTE: install via bioconda in a scanpy environment: conda create -n scanpy python=3.6 scanpy)
vii.	Bowtie 2
2.	Python libraries: math, pulp, pandas, numpy
3.	R libraries: mixdist, stringr, grDevices, RColorBrewer, pheatmap, tibble, dplyr, reshape, ggplot2, stringdist, gridExtra, KernSmooth (optional if not running on linux).
NOTE: If you don’t have these packages installed, then simply open R terminal and run: install.packages(c(“mixdist”, “stringr”, “grDevices”, “RColorBrewer”, “pheatmap”, “tibble”, “dplyr”, “reshape”, “ggplot2”,”stringdist”, “KernSmooth”,”gridExtra”))

Pipeline running procedures:
1.	Have the gene expression and AbSeq reference file ready and uploaded on Seven Bridges.
2.	Use the gene expression + AbSeq + SampleTag fastq files to run the BD Rhapsody pipeline on Seven Bridges; Wait till it finishes.
3.	Have the tetramer barcode fastq files, TCRa/b fastq files, peptide reference file and Seven Bridges DBEC csv file output and tsne coordinate csv file uploaded onto local cluster.
4.	Copy the all run scripts to the working directory:
5.	Run integrated tetramer and TCR pipelines (change file names accordingly):
bash run.sh Rhapsody.csv tetramer_R1.fastq.gz tetramer_R2.fastq.gz TCRa_R1.fastq.gz TCRa_R2.fastq.gz TCRb_R1.fastq.gz TCRb_R2.fastq.gz sampletag_R1.fastq.gz sampletag_R2.fastq.gz tSNE_coordinate.csv > log &

Explanation of output files:
Current working directory
|---- tmp (contains preliminary fastq processing files
|---- fq (contains final truecell and trimmed fastq/fasta files)
|---- tet (contains tetramer alignment and sam/bam files used for counting)
|---- out
	|---- obs.csv (the cluster information for each cell)
	|---- obsm.csv (the umap coordinate for each cell)
	|---- var.csv (number of cells detected for each genes)
	|---- marker_gene.csv (wilcox test for differential genes among clusters)
	|---- tetramer_dbec.csv (final UMI count table for tetramer tags)
	|---- tcra_final.csv (final TCRa table)
	|---- tcrb_final.csv (final TCRb table)
	|---- sampletag_umi.csv (sampletag UMI count tabe)
	|---- sampletag.csv (sampletag assignment table for each cell)
	|---- tetramer_analysis.RData (Rdata file for further customized anlaysis)
|---- figures (contain multiple QC figures)
|---- tcra (intermediate files and scripts for running TCRa data)
|---- tcrb (intermediate files and scripts for running TCRb data)
|---- truecell.txt (converted cell barcode nucleotides)
|---- tet_statistics.csv (tetramer sequencing statistics)
|---- TCRa/b_statistics.csv (TCRa/b sequencing statistics)
|---- sampletag_statistics.csv (Sampletag sequencing statistics)
|---- tetramer_dbec.csv (tetramer barcode count file)
|---- error.log (this will only output when there is error when during tetramer alignment, due to out of memory)
