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
                            os.system(
                            "metaphlan {forward_f}, {reverse_r} --input_type fastq --bowtie2db metaphlan_database --index mpa_v30_CHOCOPhlAn_201901 --bowtie2out metaphlan_output_txt/{id}.bowtie2.bz2 --tax_lev s -o metaphlan_output_txt/{id}.txt"
                            )
    else:
        id_SE = []
        for id_SE_filename in filenames:
            id_SE.append(id_SE_filename.split(".")[0])
        for id_SE_filename in filenames:
            for id in id_SE:
                if id in id_SE_filename:
                    os.system(
                    "metaphlan {id_SE_filename} --input_type fastq --bowtie2db metaphlan_database --index mpa_v30_CHOCOPhlAn_201901 --bowtie2out metaphlan_output_txt/{id}.bowtie2.bz2 --tax_lev s -o metaphlan_output_txt/{id}.txt"
                    )
    ''')
"""
}

process metaphlan_yml {
    publishDir "${params.input_directory}", mode: 'copy'

    output:
    path "metaphlan_3.0.14.yml"

    script:
    """
    #!/usr/bin/env python3
with open("metaphlan_3.0.14.yml", "w") as f:
    f.write(
        '''
name: metaphlan
channels:
  - bioconda
  - conda-forge
  - pkgs/main
dependencies:
  - _openmp_mutex=4.5=20_gnu
  - bcbio-gff=0.7.1=pyhdfd78af_3
  - biom-format=2.1.12=py37hda87dfa_1
  - biopython=1.79=py37h540881e_2
  - blast=2.16.0=hc155240_3
  - bottleneck=1.3.5=py37hda87dfa_0
  - bowtie2=2.5.4=h7071971_4
  - brotli=1.2.0=hed03a55_1
  - brotli-bin=1.2.0=hb03c661_1
  - brotli-python=1.0.9=py37hd23a5d3_7
  - bx-python=0.9.0=py37hfa8f21c_1
  - bzip2=1.0.8=hda65f42_9
  - c-ares=1.34.8=hb03c661_0
  - ca-certificates=2026.7.22=hbd8a1cb_0
  - cached-property=1.5.2=hd8ed1ab_2
  - cached_property=1.5.2=pyha770c72_1
  - capnproto=1.0.2=h2b92303_1
  - certifi=2024.8.30=pyhd8ed1ab_0
  - charset-normalizer=3.4.0=pyhd8ed1ab_0
  - click=8.1.3=py37h89c1867_0
  - cmseq=1.0.4=pyhb7b1952_0
  - curl=8.8.0=he654da7_1
  - cycler=0.11.0=pyhd8ed1ab_0
  - dendropy=5.0.11=pyhdfd78af_0
  - diamond=2.1.10=h43eeafb_2
  - entrez-direct=22.4=he881be0_0
  - fasttree=2.2.0=h7b50bb2_1
  - fonttools=4.38.0=py37h540881e_0
  - freetype=2.12.1=h267a509_2
  - future=0.18.2=py37h89c1867_5
  - gawk=5.4.1=h0a3468a_0
  - gmp=6.3.0=hac33072_2
  - gsl=2.7=he838d99_0
  - h5py=3.7.0=nompi_py37hf1ce037_101
  - hdf5=1.12.2=nompi_h4df4325_101
  - htslib=1.21=h5efdd21_0
  - icu=78.3=h54a6638_2
  - idna=3.10=pyhd8ed1ab_0
  - importlib-metadata=4.11.4=py37h89c1867_0
  - iqtree=3.1.3=h8471819_0
  - jpeg=9e=h0b41bf4_3
  - keyutils=1.6.3=hb9d3cd8_0
  - kiwisolver=1.4.4=py37h7cecad7_0
  - krb5=1.21.3=h659f571_0
  - lcms2=2.14=h6ed2654_0
  - ld_impl_linux-64=2.46.1=bootstrap_ha15bf96_2
  - lerc=4.2.0=hdb68285_0
  - libaec=1.1.5=h088129d_0
  - libasprintf=0.25.1=h3f43e3d_1
  - libblas=3.9.0=20_linux64_openblas
  - libbrotlicommon=1.2.0=hb03c661_1
  - libbrotlidec=1.2.0=hb03c661_1
  - libbrotlienc=1.2.0=hb03c661_1
  - libcblas=3.9.0=20_linux64_openblas
  - libcurl=8.8.0=hca28451_1
  - libdeflate=1.25=h17f619e_0
  - libedit=3.1.20250104=pl5321h7949ede_0
  - libev=4.33=hd590300_2
  - libffi=3.4.6=h2dba641_1
  - libgcc=16.1.0=ha9f2e26_1
  - libgcc-ng=16.1.0=h69a702a_1
  - libgettextpo=0.25.1=h3f43e3d_1
  - libgfortran=16.1.0=h69a702a_1
  - libgfortran-ng=16.1.0=h69a702a_1
  - libgfortran5=16.1.0=h79bb938_1
  - libgomp=16.1.0=he0feb66_1
  - libiconv=1.18=h3b78370_2
  - libidn2=2.3.8=hfac485b_1
  - liblapack=3.9.0=20_linux64_openblas
  - libnghttp2=1.58.0=h47da74e_1
  - libnsl=2.0.1=hb9d3cd8_1
  - libopenblas=0.3.25=pthreads_h413a1c8_0
  - libpng=1.6.43=h2797004_0
  - libsqlite=3.46.0=hde9e2c9_0
  - libssh2=1.11.0=h0841786_0
  - libstdcxx=16.1.0=h934c35e_1
  - libstdcxx-ng=16.1.0=hdf11a46_1
  - libtiff=4.4.0=h82bc61c_5
  - libunistring=0.9.10=h7f98852_0
  - libwebp-base=1.6.0=hd42ef1d_0
  - libxcb=1.13=h7f98852_1004
  - libxcrypt=4.4.36=hd590300_1
  - libzlib=1.2.13=h4ab18f5_6
  - lzo=2.10=h280c20c_1002
  - mafft=7.525=h031d066_1
  - mash=2.3=hc74b729_7
  - matplotlib-base=3.5.3=py37hf395dca_2
  - metaphlan=3.0.14=pyhb7b1952_0
  - mpfr=4.2.2=he0a73b1_0
  - munkres=1.0.7=py_1
  - muscle=5.3=h9948957_3
  - ncbi-vdb=3.4.1=hd63eeec_0
  - ncurses=6.6=hdb14827_0
  - nomkl=1.0=h5ca1d4c_0
  - numexpr=2.8.3=py37h85a3170_100
  - numpy=1.21.6=py37h976b520_0
  - openjpeg=2.5.0=h7d73246_1
  - openssl=3.6.3=h35e630c_0
  - packaging=23.2=pyhd8ed1ab_0
  - pandas=1.3.5=py37h8c16a72_0
  - patsy=0.5.6=pyhd8ed1ab_0
  - perl=5.32.1=7_hd590300_perl5
  - perl-archive-tar=3.12=pl5321hdfd78af_0
  - perl-carp=1.38=pl5321hdfd78af_4
  - perl-common-sense=3.75=pl5321hdfd78af_0
  - perl-compress-raw-bzip2=2.214=pl5321hda65f42_0
  - perl-compress-raw-zlib=2.214=pl5321h4dac143_0
  - perl-encode=3.19=pl5321hec16e2b_1
  - perl-exporter=5.72=pl5321hdfd78af_2
  - perl-exporter-tiny=1.002002=pl5321hdfd78af_0
  - perl-extutils-makemaker=7.70=pl5321hd8ed1ab_0
  - perl-io-compress=2.216=pl5321h503566f_0
  - perl-io-zlib=1.15=pl5321hdfd78af_1
  - perl-json=4.11=pl5321hdfd78af_0
  - perl-json-xs=4.04=pl5321h9948957_0
  - perl-list-moreutils=0.430=pl5321hdfd78af_0
  - perl-list-moreutils-xs=0.430=pl5321h7b50bb2_5
  - perl-parent=0.236=pl5321hdfd78af_2
  - perl-pathtools=3.75=pl5321hec16e2b_3
  - perl-scalar-list-utils=1.62=pl5321hec16e2b_1
  - perl-types-serialiser=1.01=pl5321hdfd78af_0
  - phylophlan=3.1.1=pyhdfd78af_0
  - pillow=9.2.0=py37h850a105_2
  - pip=24.0=pyhd8ed1ab_0
  - pthread-stubs=0.4=hb9d3cd8_1002
  - pyparsing=3.1.4=pyhd8ed1ab_0
  - pysam=0.16.0=py37ha9a96c6_0
  - pysocks=1.7.1=py37h89c1867_5
  - python=3.7.12=hf930737_100_cpython
  - python-dateutil=2.9.0=pyhd8ed1ab_0
  - python-lzo=1.14=py37hd06c622_1
  - python_abi=3.7=4_cp37m
  - pytz=2024.2=pyhd8ed1ab_0
  - raxml=8.2.13=h7b50bb2_3
  - readline=8.3=h853b02a_0
  - requests=2.32.2=pyhd8ed1ab_0
  - rpsbproc=0.5.0=h6a68c12_1
  - samtools=1.21=h50ea8bc_0
  - scipy=1.7.3=py37hf2a6cf1_0
  - seaborn=0.12.2=hd8ed1ab_0
  - seaborn-base=0.12.2=pyhd8ed1ab_0
  - setuptools=69.0.3=pyhd8ed1ab_0
  - six=1.16.0=pyh6c4a22f_0
  - sqlite=3.46.0=h6d4b2fc_0
  - statsmodels=0.13.2=py37hda87dfa_0
  - tk=8.6.13=noxft_h4845f30_101
  - trimal=1.5.1=h9948957_0
  - typing-extensions=4.7.1=hd8ed1ab_0
  - typing_extensions=4.7.1=pyha770c72_0
  - unicodedata2=14.0.0=py37h540881e_1
  - urllib3=2.2.1=pyhd8ed1ab_0
  - wget=1.21.4=hda4d442_0
  - wheel=0.42.0=pyhd8ed1ab_0
  - xorg-libxau=1.0.12=hb03c661_1
  - xorg-libxdmcp=1.1.5=hb03c661_1
  - xz=5.2.6=h166bdaf_0
  - zipp=3.15.0=pyhd8ed1ab_0
  - zlib=1.2.13=h4ab18f5_6
  - zstd=1.5.6=ha6fb4c9_0
''')
    """
}

process metaphlan {

    publishDir "${params.input_directory}", mode: 'copy'

    input:
    path x

    output:
    path "metaphlan_output_txt"
    path "metaphlan_merged_table.tsv"
    
    script:
    '''
    conda env create -f metaphlan_3.0.14.yml
    mkdir metaphlan_output_txt
    mkdir metaphlan_database
    mamba run -n metaphlan metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db metaphlan_database
    mamba run -n metaphlan python3 ${params.input_directory}metaphlan_command.py
    mamba run -n metaphlan metaphlan merge_metaphlan_tables.py metaphlan_merged_table.tsv > metaphlan_merged_table.tsv 
    '''
}


workflow {
    req_fastqs = Channel.fromPath("${params.input_directory}/*")
    fastqc(req_fastqs)
    metaphlan_script()
    metaphlan_yml()
    metaphlan(metaphlan_yml.out)
}
