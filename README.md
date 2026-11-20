The given repo contains two separate nextflow based pipelines for downstream analysis of 16S and Whole Genome Sequencing (WGS) reads using QIIME2 and Metaphlan3 respectively. 
In both of the pipelines, the analysis is precceded by FASTQC analysis. 
The main selling point of these pipelines are the fact that they require minimal human intervention, with only the required input being the path to the directory containing the FASTQ files. 

Please ensure that your input directory contains only FASTQ files. If any other file types are detected, the pipeline will not run.

In the 16S pipeline, the following steps constitute the pipeline:
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
For more important information regaridng the 16S pipeline, please visit:
https://github.com/Rohithposts/Nextflow_pipelines_16S_WGS/wiki/16S-Pipeline


Similarly, the WGS nextflow pipeline accepts the path to the WGS fastq files as input and runs FASTQC and Metaphlan3, storing the final results in the input directory.
```bash
nextflow run nextflow_WGS.nf --input_directory /Path/to/Your/WGS/Fastq/Files/ -with-conda
```
For more important information regaridng the WGS pipeline, please visit:
https://github.com/Rohithposts/Nextflow_pipelines_16S_WGS/wiki/WGS-Pipeline

Importantly, when running these pipelines, please ensure that your input directory contains only FASTQ files. If any other files are also detected, the pipeline doesnt run.
Please ensure that you have an active internet connection when running the pipeline, as the pipeline downloads mutliplte conda based tools.
