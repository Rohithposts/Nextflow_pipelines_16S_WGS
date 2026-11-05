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

process metaphlan_script {

    publishDir "${params.input_directory}", mode: 'copy'

    output:
    path "metaphlan_command.py"

    script:
    """
    #!/usr/bin/env python3
    with open("metaphlan_command.py", "w") as f:
        f.write(
    '''
    #!/usr/bin/env python3
    import os
    from pathlib import Path
    input_dir = Path("${params.input_directory}")
    filenames = []
    for file in input_dir.iterdir():
        filenames.append(file.name)
    if "_" in filenames[0]:
        ids_PE = []
        for x in filenames:
            ids_PE.append(x.split("_")[0])
        PE_f = ["_1", "_R1"]
        PE_r = ["_2", "_R2"]
        for f in PE_f:
            for r in PE_r:
                for id in ids_PE:
                    for file in filenames:
                        if id + f in file:
                            forward_f = file
                        if id + r in file:
                            reverse_r = file
                            os.system("
                            metaphlan {forward_f}, {reverse_r} \
                            --input type fastq \
                            --bowtie2db metaphlan_database \
                            --index mpa_v30_CHOCOPhlAn_201901 \
                            --bowtie2out metaphlan_output_txt/{id}.bowtie2.bz2 \
                            --tax_lev s \
                            -o metaphlan_output_txt/{id}.txt; done
                            ")
    else:
        id_SE = []
        for id_SE_filename in filenames:
            id_SE.append(id_SE_filename.split(".")[0])
        for id_SE_filename in filenames:
            for id in id_SE:
                if id in id_SE_filename:
                    os.system("
                    metaphlan {id_SE_filename} \
                    --input type fastq \
                    --bowtie2db metaphlan_database \
                    --index mpa_v30_CHOCOPhlAn_201901 \
                    --bowtie2out metaphlan_output_txt/{id}.bowtie2.bz2 \
                    --tax_lev s \
                    -o metaphlan_output_txt/{id}.txt
                    ")
    ''')
"""
}

process metaphlan {

    conda 'conda-forge::python=3.10 bioconda::metaphlan=3.0'

    publishDir "${params.input_directory}", mode: 'copy'


    input:
    path reads_2

    output:
    path "metaphlan_output_txt"
    path "metaphlan_merged_table.tsv"
    script:
    """
    mkdir "metaphlan_output_txt"
    mkdir "metaphlan_database"
    metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db metaphlan_database
    python3 ${params.input_directory}metaphlan_command.py
    metaphlan merge_metaphlan_tables.py metaphlan_merged_table.tsv > metaphlan_merged_table.tsv 
    """
}


workflow {
    req_fastqs = Channel.fromPath("${params.input_directory}/*")
    fastqc(req_fastqs)
    metaphlan(req_fastqs)
}
