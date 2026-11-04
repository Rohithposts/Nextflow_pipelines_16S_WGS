nextflow.enable.dsl=2
params.input_directory = null



process fastqc {
    conda 'bioconda::fastqc=0.12.1'

    publishDir "${params.input_directory}", mode: 'copy'

    input:
    val reads

    output:
    path "fastqc_results/*"

    script:
    """
    mkdir -p fastqc_results
    fastqc ${reads} -o fastqc_results

    """
}

process QIIME_manifest {
    publishDir "${params.input_directory}", mode: 'copy'

    input:
    val reads_ch

    output:
    path "manifest_*.tsv"

    script:
    """
    #!/usr/bin/env python3
from pathlib import Path
import os
cur_dir = Path("${reads_ch}")
file_names = []
for file in cur_dir.iterdir():
    file_names.append(file.name)
    print(file.name)
file1 = file_names[0]
print(file1)
if "_" in file1:
    ids = []
    fwd = []
    rvrse = []
    for file in cur_dir.iterdir():
        ids.append(file.name.split("_")[0])
        if "_1" in file.name or "_R1" in file.name:
            fwd.append(os.path.abspath(file))
        if "_2" in file.name or "_R2" in file.name:
            rvrse.append(os.path.abspath(file))
        with open("manifest_PE.tsv", "w") as f:
            f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\\n")
            for x in set(ids):
                for y in fwd:
                    for z in rvrse:
                        if x in y:
                            if x in z:
                                f.write(f"{x}\\t{y}\\t{z}\\n")
else:
    ids = []
    filepaths = []
    for file in cur_dir.iterdir():
        ids.append(file.name.split(".")[0])
        filepaths.append(os.abspath.path(file))
    with open("manifest_SE.tsv", "w") as f:
        f.write("sample-id\tabsolute-filepath\\n")
    for x in set(ids):
        for y in filepaths:
            if x in y:
                f.write(f"{x}\\t{y}\\n")

"""
}


process primer_detection {

    publishDir "${params.input_directory}", mode: 'copy'
    input:
    path fastq_file
    path manifest
    output:
    path "primer_detected.txt"
    path "primer_detected_summary.txt"
    script:
    """
    #!/usr/bin/env python3
import re
import gzip
import os
from pathlib import Path
from collections import Counter
dir = os.path.abspath("${fastq_file}")
manifest = Path("${manifest}")
file_names_collection = []
for file in Path("${params.input_directory}").iterdir():
    if not file.name == "fastqc_results":
        file_names_collection.append(file.name)
if manifest.name == "manifest_SE.tsv":
    fastqfile_names = []
    for filename in Path("${params.input_directory}").iterdir():
        if not filename.name == "fastqc_results":
            if not filename.name == "manifest_SE.tsv":
                fastqfile_names.append(filename) 
    forward_primers = ["AGAGTTTGAT[CT][AC]TGGCTCAG", "GGGTTCGATTCTGGCTCAG", "GACTCCTACGGGAGGC[AT]GCAG", "ACTCCTACGGGAGGC[AT]GCAG", "CCTACGGG[ATGC]GGC[AT]GCAG", "GTG[CT]CAGC[AC]GCCGCGGTAA", "[CT]CAGC[AC]GCCGCGGTAATACG", "[AG]GGATTAGATACCC", "GGATTAGATACCCTGGTAGTCC", "GAATTGACGGGGGCCCGCACAAG", "AACT[CT]AAA[GT]GAATTG[AG]CGG", "CAACGAGCGCAACCCT"]
    reverse_primers = ["CTGC[AT]GCCTCCCGTAGGAGTC", "GC[AT]GCCTCCCGTAGGAGT", "CGTATTACCGCGGCTGCTGG", "ATTACCGCGGCTGCTGG", "GACTAC[ACT][ACG]GGGT[AT]TCTAATCC", "TAC[ACT]AGGGT[AT]TCTAATCC", "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT", "GAATTAAACCACATGCTC", "CCG[CT]CAATT[CT][AC]TTT[AG]AGTTT", "CGAC[AG][AG]CCATGCA[ATGC]CACCT", "CGGTGTGTACAAGGCCCGGGAACG", "TACGG[CT]TACCTTGTTACGACTT"]
    forwardprimers_detected = Counter()
    reverseprimers_detected = Counter()
    def temp_iupac_cnvtr(primer):
        iupac_ambiguity_function = {
        "[AG]": "R",
        "[CT]": "Y",
        "[CG]": "S",
        "[AT]": "W",
        "[GT]": "K",
        "[AC]": "M",
        "[CGT]": "B",
        "[AGT]": "D",
        "[ACT]": "H",
        "[ACG]": "V",
        "[ACGT]": "N"
        }
        output_iupac = ""
        for x, y in iupac_ambiguity_function.items():
            primer = primer.replace(x,y)
        output_iupac = primer
        return output_iupac
    def reverse_compliment(primer):
        rev_comp = []
        revs = {
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G',
        'H':'D',
        'D':'H',
        'Y':'R',
        'R':'Y',
        'S':'S',
        'W':'W',
        'K':'M',
        'M':'K',
        'B':'V',
        'V':'B',
        'D':'H',
        'H':'D',
        'N':'N'
        }
        for letter in reversed(primer):
            rev_comp.append(revs[letter])
        return ''.join(rev_comp)
    def final_rev_search_input(primer_3):
        primer_A = reverse_compliment(temp_iupac_cnvtr(primer_3))
        iupac_ambiguity_function_final = {
        "R": "[AG]",
        "Y": "[CT]",
        "S": "[CG]",
        "W": "[AT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
        "N": "[ACGT]"
        }
        output_iupac_A = ""
        for x, y in iupac_ambiguity_function_final.items():
            primer_A = primer_A.replace(x,y)
        output_iupac_A = primer_A
        return output_iupac_A
    for fastq_file_PD in fastqfile_names:
        if fastq_file_PD.name.endswith("gz"):
            with gzip.open(fastq_file_PD, "rt") as f:  
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        for primer in forward_primers:
                            if re.search(primer,line.strip()[:50]):
                                forwardprimers_detected[primer] += 1
                        for primer in reverse_primers:
                            if re.search(final_rev_search_input(primer), line.strip()[:50]):
                                reverseprimers_detected[primer] += 1

        else:
            with open(fastq_file_PD) as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        for primer in forward_primers:
                            if re.search(primer,line.strip()[:50]):
                                forwardprimers_detected[primer] += 1
                        for primer in reverse_primers:
                            if re.search(final_rev_search_input(primer), line.strip()[:50]):
                                reverseprimers_detected[primer] += 1
                                
    iupac_ambiguity = {
    "[AG]": "R",
    "[CT]": "Y",
    "[CG]": "S",
    "[AT]": "W",
    "[GT]": "K",
    "[AC]": "M",
    "[CGT]": "B",
    "[AGT]": "D",
    "[ACT]": "H",
    "[ACG]": "V",
    "[ACGT]": "N"
    }
    new_forwardprimers_detected = {}
    new_reverseprimers_detected = {}
    for primer_f, detected_count_f in forwardprimers_detected.items():
        for ambigous, unambigous in iupac_ambiguity.items():
            primer_f = primer_f.replace(ambigous, unambigous)
        new_forwardprimers_detected[primer_f] = detected_count_f
    for primer_r, detected_count_r in reverseprimers_detected.items():
        for ambigous, unambigous in iupac_ambiguity.items():
            primer_r = primer_r.replace(ambigous, unambigous)
        new_reverseprimers_detected[primer_r] = detected_count_r
    fprimer = ""
    rprimer = ""
    fprimer_counts = []
    rprimer_counts = []
    for value, count in new_forwardprimers_detected.items():
        fprimer_counts.append(count)
    for value, count in new_reverseprimers_detected.items():
        rprimer_counts.append(count)
    for value, count in new_forwardprimers_detected.items():
        if count == max(fprimer_counts):
            fprimer = value
    for value, count in new_reverseprimers_detected.items():
        if count == max(rprimer_counts):
            rprimer = value    
    with open("primer_detected.txt", "w") as f:
        f.write(fprimer + "\\n")
        f.write(rprimer + "\\n")
    with open("primer_detected_summary.txt", "w") as f:
        f.write("List of all primers detected:\\n")
        f.write("Forward\\n")
        for x,y in new_forwardprimers_detected.items():
            f.write(f"{x} was detected {y} times\\n")
        f.write("Reverse\\n")
        for z,w in new_reverseprimers_detected.items():
            f.write(f"{z} was detected {w} times\\n")
        f.write("List of FASTQ Files searched for primer detection:\\n")
        for searchedfile in fastqfile_names:
            f.write(searchedfile + "\\n")
elif manifest.name == "manifest_PE.tsv":
    fastq_file_collection = []
    id_s = []
    for file in Path("${params.input_directory}").iterdir():
        if not file.name == "fastqc_results":
            if not "manifest_PE" in file.name:
                fastq_file_collection.append(file.name)
                id_s.append(file.name.split("_")[0])
    reqids = set(id_s)
    files_searched = []
    primer_detected_f = []
    primer_detected_r = []
    primercounts_f = Counter()
    primercounts_r = Counter()
    for reqid in reqids: 
        forwardfile = ""
        reversefile = ""
        file_endings_forward = ["_1", "_R1"]
        file_endings_reverse = ["_2", "_R2"]
        for file in Path("${params.input_directory}").iterdir():
            if not file.name == "fastqc_results":
                if not "manifest" in file.name:
                    for forward in file_endings_forward:
                        if reqid + forward in file.name:
                            forwardfile = file
                            files_searched.append(file.name)
                    for reverse in file_endings_reverse:
                        if reqid + reverse in file.name:
                            reversefile = file
                            files_searched.append(file.name)
        forward_primers = ["AGAGTTTGAT[CT][AC]TGGCTCAG", "GGGTTCGATTCTGGCTCAG", "GACTCCTACGGGAGGC[AT]GCAG", "ACTCCTACGGGAGGC[AT]GCAG", "CCTACGGG[ATGC]GGC[AT]GCAG", "GTG[CT]CAGC[AC]GCCGCGGTAA", "[CT]CAGC[AC]GCCGCGGTAATACG", "[AG]GGATTAGATACCC", "GGATTAGATACCCTGGTAGTCC", "GAATTGACGGGGGCCCGCACAAG", "AACT[CT]AAA[GT]GAATTG[AG]CGG", "CAACGAGCGCAACCCT"]
        if forwardfile.name.endswith("gz"):
            with gzip.open(forwardfile, "rt") as f:  
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        for primer in forward_primers:
                            if re.search(primer,line.strip()[:50]):
                                primercounts_f[primer] += 1
                                
        else:
            with open(forwardfile) as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        for primer in forward_primers:
                            if re.search(primer,line.strip()[:50]):
                                primercounts_f[primer] += 1
                            
        reverse_primers = ["CTGC[AT]GCCTCCCGTAGGAGTC", "GC[AT]GCCTCCCGTAGGAGT", "CGTATTACCGCGGCTGCTGG", "ATTACCGCGGCTGCTGG", "GACTAC[ACT][ACG]GGGT[AT]TCTAATCC", "TAC[ACT]AGGGT[AT]TCTAATCC", "GGACTAC[ACT][ACG]GGGT[AT]TCTAAT", "GAATTAAACCACATGCTC", "CCG[CT]CAATT[CT][AC]TTT[AG]AGTTT", "CGAC[AG][AG]CCATGCA[ATGC]CACCT", "CGGTGTGTACAAGGCCCGGGAACG", "TACGG[CT]TACCTTGTTACGACTT"]                        
        if reversefile.name.endswith("gz"):
            with gzip.open(reversefile, "rt") as f:  
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        for primer in reverse_primers:
                            if re.search(primer,line.strip()[:50]):
                                primercounts_r[primer] += 1
                            
        else:
            with open(reversefile) as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        for primer in reverse_primers:
                            if re.search(primer,line.strip()[:50]):
                                primercounts_r[primer] += 1
    iupac_ambiguity = {
    "[AG]": "R",
    "[CT]": "Y",
    "[CG]": "S",
    "[AT]": "W",
    "[GT]": "K",
    "[AC]": "M",
    "[CGT]": "B",
    "[AGT]": "D",
    "[ACT]": "H",
    "[ACG]": "V",
    "[ACGT]": "N"
    }
    new_primercounts_f = {}
    new_primercounts_r = {}
    for primer_f, detected_count_f in primercounts_f.items():
        for ambigous, unambigous in iupac_ambiguity.items():
            primer_f = primer_f.replace(ambigous, unambigous)
        new_primercounts_f[primer_f] = detected_count_f
    for primer_r, detected_count_r in primercounts_r.items():
        for ambigous, unambigous in iupac_ambiguity.items():
            primer_r = primer_r.replace(ambigous, unambigous)
        new_primercounts_r[primer_r] = detected_count_r
    fprimer = ""
    rprimer = ""
    fprimer_counts = []
    rprimer_counts = []
    for value, count in new_primercounts_f.items():
        fprimer_counts.append(count)
    for value, count in new_primercounts_r.items():
        rprimer_counts.append(count)
    for value, count in new_primercounts_f.items():
        if count == max(fprimer_counts):
            fprimer = value
    for value, count in new_primercounts_r.items():
        if count == max(rprimer_counts):
            rprimer = value    
    with open("primer_detected.txt", "w") as f:
        f.write(fprimer + "\\n")
        f.write(rprimer + "\\n")
    with open("primer_detected_summary.txt", "w") as f:
        f.write("List of all primers detected:\\n")
        f.write("Forward\\n")
        for x,y in new_primercounts_f.items():
            f.write(f"{x} was detected {y} times\\n")
        f.write("Reverse\\n")
        for z,w in new_primercounts_r.items():
            f.write(f"{z} was detected {w} times\\n")
        f.write("List of FASTQ Files searched for primer detection:\\n")
        for searchedfile in files_searched:
            f.write(searchedfile + "\\n")
    """
}

process check_phred_quality_coding {
    publishDir "${params.input_directory}", mode: 'copy'

    input:
    path read

    output:
    path "manifest_*_Phred*.tsv"

    script:
    """
    #!/usr/bin/env python3
    import gzip
    from pathlib import Path
    import os
    import shutil
    qs = ""
    letters_ascii = []
    filenames = []
    inputdir = Path("${params.input_directory}")
    path_quality_check = Path("${params.input_directory}")
    read_absdir = path_quality_check.resolve()
    for x in path_quality_check.iterdir():
        if not x.name == "fastqc_results":
            if not x.name.endswith("tsv"):
                filenames.append(x)
    file = filenames[0]
    if file.name.endswith(".gz"):
        with gzip.open(file, "rt") as f:
            for index, line in enumerate(f):
                if index == 3:
                    for letter in line:
                        letters_ascii.append(ord(letter))
            if min(letters_ascii) < 59:
                qs = "Phred33V2"
            else:
                qs = "Phred64V2"
    else:
        with open(file) as f:
            for index, line in enumerate(f):
                if index == 3:
                    for letter in line:
                        letters_ascii.append(ord(letter))
            if min(letters_ascii) < 59:
                qs = "Phred33V2"
            else:
                qs = "Phred64V2"
    for file in inputdir.iterdir():
        if file.name.startswith("manifest_"):
            req_file = file
            req_file = file.with_name(f"{file.stem}_{qs}.tsv")
            file.rename(req_file)
            shutil.copy(req_file, '.')
"""
}


process qiime {    

    publishDir "${params.input_directory}", mode: 'copy'


    output:
        path "QIIME_Artifacts"


    script:
    """
    mkdir -p QIIME_Artifacts
    wget https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.10/qiime2/released/rachis-qiime2-linux-64-conda.yml
    conda env create -f rachis-qiime2-linux-64-conda.yml -n qiime2
    primer1=\$(head -n 1 "${params.input_directory}primer_detected.txt") 
    primer2=\$(sed -n '2p' "${params.input_directory}primer_detected.txt")
    
    if [[ -f "${params.input_directory}manifest_SE_Phred33V2.tsv" ]]; then
        mamba run -n qiime2 qiime tools import --type 'SampleData[SequencesWithQuality]' \
        --input-path ${params.input_directory}manifest_SE_Phred33V2.tsv \
        --output-path QIIME_Artifacts/samples.qza \
        --input-format SingleEndFastqManifestPhred33V2
        mamba run -n qiime2 qiime cutadapt trim-single --i-demultiplexed-sequences QIIME_Artifacts/samples.qza \
        --p-front-f "\$primer1" \
        --p-front-f "\$primer2" \
        --p-front-r "\$primer2" \
        --p-front-r "\$primer1" \
        --p-match-adapter-wildcards --p-no-discard-untrimmed --p-cores 14 --o-trimmed-sequences QIIME_Artifacts/primers_removed.qza --o-stats QIIME_Artifacts/primer-trimming-stats.qza
        mamba run -n qiime2 qiime dada2 denoise-single --i-demultiplexed-seqs QIIME_Artifacts/primers_removed.qza --p-trunc-len-f 0 --p-trunc-len-r 0 \
        --o-table QIIME_Artifacts/table.qza --o-representative-sequences QIIME_Artifacts/seqs.qza \
        --o-base-transition-stats QIIME_Artifacts/base-transition-stats.qza \
        --o-denoising-stats QIIME_Artifacts/stats.qza --p-n-threads 10
    
    fi

    if [[ -f "${params.input_directory}manifest_PE_Phred33V2.tsv" ]]; then
        mamba run -n qiime2 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${params.input_directory}manifest_PE_Phred33V2.tsv \
        --output-path QIIME_Artifacts/samples.qza \
        --input-format PairedEndFastqManifestPhred33V2
        mamba run -n qiime2 qiime cutadapt trim-paired --i-demultiplexed-sequences QIIME_Artifacts/samples.qza \
        --p-front-f "\$primer1" \
        --p-front-f "\$primer2" \
        --p-front-r "\$primer2" \
        --p-front-r "\$primer1" \
        --p-match-adapter-wildcards --p-no-discard-untrimmed --p-cores 14 --o-trimmed-sequences QIIME_Artifacts/primers_removed.qza --o-stats QIIME_Artifacts/primer-trimming-stats.qza
        mamba run -n qiime2 qiime dada2 denoise-paired --i-demultiplexed-seqs QIIME_Artifacts/primers_removed.qza --p-trunc-len-f 0 --p-trunc-len-r 0 \
        --o-table QIIME_Artifacts/table.qza --o-representative-sequences QIIME_Artifacts/seqs.qza \
        --o-base-transition-stats QIIME_Artifacts/base-transition-stats.qza \
        --o-denoising-stats QIIME_Artifacts/stats.qza --p-n-threads 10

    fi

    if [[ -f "${params.input_directory}manifest_SE_Phred64V2.tsv" ]]; then
        mamba run -n qiime2 qiime tools import --type 'SampleData[SequencesWithQuality]' \
        --input-path ${params.input_directory}manifest_SE_Phred64V2.tsv \
        --output-path QIIME_Artifacts/samples.qza \
        --input-format SingleEndFastqManifestPhred64V2
        mamba run -n qiime2 qiime cutadapt trim-single --i-demultiplexed-sequences QIIME_Artifacts/samples.qza \
        --p-front-f "\$primer1" \
        --p-front-f "\$primer2" \
        --p-front-r "\$primer2" \
        --p-front-r "\$primer1" \
        --p-match-adapter-wildcards --p-no-discard-untrimmed --p-cores 14 --o-trimmed-sequences QIIME_Artifacts/primers_removed.qza --o-stats QIIME_Artifacts/primer-trimming-stats.qza
        mamba run -n qiime2 qiime dada2 denoise-single --i-demultiplexed-seqs QIIME_Artifacts/primers_removed.qza --p-trunc-len-f 0 --p-trunc-len-r 0 \
        --o-table QIIME_Artifacts/table.qza --o-representative-sequences QIIME_Artifacts/seqs.qza \
        --o-base-transition-stats QIIME_Artifacts/base-transition-stats.qza \
        --o-denoising-stats QIIME_Artifacts/stats.qza --p-n-threads 10

    fi

    if [[ -f "${params.input_directory}manifest_PE_Phred64V2.tsv" ]]; then
        mamba run -n qiime2 qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${params.input_directory}manifest_PE_Phred64V2.tsv \
        --output-path QIIME_Artifacts/samples.qza \
        --input-format PairedEndFastqManifestPhred64V2
        mamba run -n qiime2 qiime cutadapt trim-paired --i-demultiplexed-sequences QIIME_Artifacts/samples.qza \
        --p-front-f "\$primer1" \
        --p-front-f "\$primer2" \
        --p-front-r "\$primer2" \
        --p-front-r "\$primer1" \
        --p-match-adapter-wildcards --p-no-discard-untrimmed --p-cores 14 --o-trimmed-sequences QIIME_Artifacts/primers_removed.qza --o-stats QIIME_Artifacts/primer-trimming-stats.qza
        mamba run -n qiime2 qiime dada2 denoise-paired --i-demultiplexed-seqs QIIME_Artifacts/primers_removed.qza --p-trunc-len-f 0 --p-trunc-len-r 0 \
        --o-table QIIME_Artifacts/table.qza --o-representative-sequences QIIME_Artifacts/seqs.qza \
        --o-base-transition-stats QIIME_Artifacts/base-transition-stats.qza \
        --o-denoising-stats QIIME_Artifacts/stats.qza --p-n-threads 10

    fi
        
    wget https://www.arb-silva.de/fileadmin/silva_databases/current/QIIME2/2026.7/SSU/full-length/uniform/SILVA_144_SSURef_NR99_full-length_prokaryotes_rank-propagation_uniq_genus_uniform_classifier.qza
    mamba run -n qiime2 qiime feature-classifier classify-sklearn --i-reads QIIME_Artifacts/seqs.qza --i-classifier SILVA_144_SSURef_NR99_full-length_prokaryotes_rank-propagation_uniq_genus_uniform_classifier.qza --p-n-jobs 3 --o-classification QIIME_Artifacts/taxonomy_classification.qza
    mamba run -n qiime2 qiime phylogeny align-to-tree-mafft-fasttree --i-sequences QIIME_Artifacts/seqs.qza --o-alignment QIIME_Artifacts/alignedseqs.qza --o-masked-alignment QIIME_Artifacts/masked-alignment.qza --o-tree QIIME_Artifacts/unrootedtree.qza --o-rooted-tree QIIME_Artifacts/rooted-tree.qza --p-n-threads 12
    mamba run -n qiime2 qiime tools export --input-path QIIME_Artifacts/table.qza --output-path QIIME_Artifacts/exported_table
    mamba run -n qiime2 biom convert -i QIIME_Artifacts/exported_table/feature-table.biom -o QIIME_Artifacts/abundance_table.tsv --to-tsv
    conda env remove -n qiime2 -y
"""
}



workflow {
    req_fastqs = Channel.fromPath("${params.input_directory}/*")
    primer_testing = req_fastqs.first()
    fastqc(req_fastqs)
    //multiqc(fastqc.out)
    QIIME_manifest(Channel.value(file(params.input_directory)))
    primer_detection(Channel.value(file(params.input_directory)), QIIME_manifest.out)
    check_phred_quality_coding(Channel.value(file(params.input_directory)))
    qiime()
}


