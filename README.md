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


Some important points to keep in mind when using the 16S pipeline:

1)The pipeline feeds the most frequently detected primer(s) to qiime cutadapt, which, in some cases, might not be the primer pair in mind. It also created a txt file listing how many detections for all primers were observed. The python re.search() function was used for this process.

2)Since the primers are already trimmed off usign qiime cutadapt, all the input numeric values for qiime dada2 are 0.

3)The SILVA Database used for taxonomy classification is the full-length,uniform SILVA_144_SSURef_NR99 classifier from the official SILVA Database.

Some important points to keep in mind when using the WGS pipeline:

1)Metaphlan3, a slightly older version is used in this pipeline. This is done as it requires a relatively smaller sized database, as opposed to a 40 GB database for metaphlan4. This is more practically feasible for more devices of different storage spaces. For taxonomic profiling, the mpa_v30_CHOCOPhlAn_201901 metaphlan database is used.

2)In the given pipeline, for all the FASTQ Files, Metaphlan analysis is done at the species level.





Contact: sreesairohith@gmail.com

