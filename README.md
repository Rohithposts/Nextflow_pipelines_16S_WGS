The given repo contains two separate nextflow based pipelines for downstream analysis of 16S and Whole Genome Sequencing (WGS) reads using QIIME2 and Metaphlan3 respectively. 
In both of the pipelines, the analysis is precceded by FASTQC analysis. 
The main selling point of these pipelines are the fact that they require minimal human intervention, with only the required input being the path to the directory containing the FASTQ files. In the 16S pipeline, the following steps constitute the pipeline:
```bash
FASTQC
QIIME Manifest Generation and checking if reads are Single End or Paired End
Primer detection
Fastq file phred encoding check
Primer trimming using QIIME cutadapt
DADA2 Denoising using qiime dada2
```
The pipeline, in itself, handles all conda related dependencies and also downlaods the required SILVA classifier; all the final files such as the abundance table, QIIME artifacts, FASTQC results etc. are stored in the same input direcotory.
The pipeline can be run with the following command:
```bash
nextflow run nextflow_16S.nf --input_directory /Path/to/Your/16S/Fastq/Files/ -with-conda
```
Similarly, the WGS nextflow pipeline accepts the path to the WGS fastq files as input and runs FASTQC and Metaphlan3, storing the final results in the input directory.
```bash
nextflow run nextflow_WGS.nf --input_directory /Path/to/Your/16S/Fastq/Files/ -with-conda
```
