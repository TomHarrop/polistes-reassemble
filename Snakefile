#!/usr/bin/env python3

import multiprocessing
import pandas
import pathlib2


#############
# FUNCTIONS #
#############

def meraculous_input_resolver(wildcards):
    # which readset are we on
    my_readset = wildcards.read_set
    if my_readset == 'norm':
        std_path = 'output/020_norm/{library}-norm.fq.gz'
    elif my_readset == 'raw':
        std_path = 'output/010_trim-decon/{library}.fq.gz'
    # mp libs are always the same
    mp_path = 'output/010_trim-decon/{library}.fq.gz'
    lib_dict = dict.fromkeys(all_libs)
    # generate paths from keys
    for key in lib_dict:
        if 'standard' in key:
            lib_dict[key] = std_path.format(library=key)
        if 'mate' in key:
            lib_dict[key] = mp_path.format(library=key)
    return(lib_dict)


def resolve_path(x):
    mypath = pathlib2.Path(x).resolve()
    return str(mypath)


def trim_resolver(wildcards):
    my_srr = name_to_srr[wildcards.library]
    return({
        'r1': f'output/fastq_repaired/{my_srr}_1.fastq.gz',
        'r2': f'output/fastq_repaired/{my_srr}_2.fastq.gz'})


###########
# GLOBALS #
###########

run_info_file = 'data/SraRunInfo.csv'
meraculous_config_file = 'src/meraculous_config.txt'

sra_container = 'shub://TomHarrop/singularity-containers:sra_2.9.2'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
mer_container = 'shub://TomHarrop/singularity-containers:meraculous_2.2.6'


########
# MAIN #
########

run_info = pandas.read_csv(run_info_file)

# split the run info into name_to_url dict
col_to_exp = run_info.to_dict()['Run']
col_to_url = run_info.to_dict()['download_path']
col_to_name = run_info.to_dict()['LibraryName']
col_to_insert = run_info.to_dict()['InsertSize']
col_to_sd = run_info.to_dict()['InsertDev']

# make the required dictionaries
name_to_url = {}
for key in col_to_exp:
    name_to_url[col_to_exp[key]] = col_to_url[key]

name_to_srr = {}
for key in col_to_name:
    my_name = col_to_name[key].replace(" ", "_").replace("_library", "")
    name_to_srr[my_name] = col_to_exp[key]

# get a list of all names
all_samples = sorted(set(name_to_url.keys()))
all_libs = sorted(set(name_to_srr.keys()))

# work out which ones are mate pairs
mate_pair_libs = {k: name_to_srr[k] for k in name_to_srr.keys()
                  if 'mate' in k}

standard_libs = {k: name_to_srr[k] for k in name_to_srr.keys()
                 if 'mate' not in k}


# read the meraculous config
with open(meraculous_config_file, 'rt') as f:
    meraculous_config_string = ''.join(f.readlines())

# format using
# meraculous_config_string.format(**lib_dict)

# fix the headers (do during trim-decon)
# reformat.sh \
# in=SRR1393722_1.fastq.gz \
# in2=SRR1393722_2.fastq.gz \
# out=stdout.fastq \
#     | head -n 16 \
#     | sed -e 's/^@\S*\s/@/g' \
#     | reformat.sh \
#         in=stdin.fastq \
#         int=t \
#         out=stdout.fastq \
#         addcolon=t \
#         trimreaddescription=t \
#     | sed -e '/^@\S\+:/s/$/N:0:NNNNNN/g'


 
#########
# RULES #
#########

rule target:
    input:
        expand(('output/040_meraculous/'
                '{read_set}_k{k}_diplo{diplo}/'
                'meraculous_final_results/final.scaffolds.fa'),
               read_set=['norm', 'raw'],
               k=[35], # , 61, 91], # have to disable MP libraries for k>35
               diplo=['0', '1'])

# assembly rule
rule meraculous:
    input:
        unpack(meraculous_input_resolver),
        config = ('output/040_meraculous/'
                  '{read_set}_k{k}_diplo{diplo}/config.txt')
    output:
        contigs = ('output/040_meraculous/'
                   '{read_set}_k{k}_diplo{diplo}/'
                   'meraculous_final_results/final.scaffolds.fa'),
    params:
        outdir = 'output/040_meraculous/{read_set}_k{k}_diplo{diplo}/',
    threads:
        multiprocessing.cpu_count()
    log:
        'output/logs/040_meraculous/{read_set}_k{k}_diplo{diplo}.log'
    singularity:
        mer_container
    shell:
        'run_meraculous.sh '
        '-dir {params.outdir} '
        '-config {input.config} '
        '-cleanup_level 2 '
        '&> {log}'

rule meraculous_config:
    input:
        unpack(meraculous_input_resolver)
    output:
        config = ('output/040_meraculous/'
                  '{read_set}_k{k}_diplo{diplo}/config.txt')
    params:
        dmin = '0'
    threads:
        multiprocessing.cpu_count()
    run:
        my_lib_dict = meraculous_input_resolver(wildcards)
        # fill in full paths
        for key in my_lib_dict:
            my_lib_dict[key] = resolve_path(my_lib_dict[key])
        my_lib_dict['k'] = wildcards.k
        my_lib_dict['diplo'] = wildcards.diplo
        my_lib_dict['dmin'] = params.dmin
        my_lib_dict['threads'] = threads
        my_conf = meraculous_config_string.format(**my_lib_dict)
        with open(output.config, 'wt') as f:
            f.write(my_conf)

# split nextera for the mate pair libs
rule split_nextera:
    input:
        fq = 'output/010_trim-decon/{library}.fq.gz',
    output:
        fq = 'output/030_split/{library}-lmp.fq.gz',
        frag = 'output/030_split/{library}-frag.fq.gz',
        sing = 'output/030_split/{library}-sing.fq.gz',
        unk = 'output/030_split/{library}-unk.fq.gz'
    log:
        'output/logs/030_split/{library}-split.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbduk_container
    shell:
        'splitnextera.sh '
        'mask '
        'interleaved '
        'ziplevel=9 '
        'in={input.fq} '
        'out={output.fq} '
        'outf={output.frag} '
        'outu={output.unk} '
        'outs={output.sing} '
        '&> {log}'

# normalise the standard libs
rule bbnorm:
    input:
        fq = 'output/010_trim-decon/{library}.fq.gz',
    output:
        fq_norm = 'output/020_norm/{library}-norm.fq.gz',
        fq_toss = 'output/020_norm/{library}-toss.fq.gz',
        hist = 'output/020_norm/{library}-hist.txt',
        hist_out = 'output/020_norm/{library}-hist_out.txt',
        peaks = 'output/020_norm/{library}-peaks.txt'
    log:
        norm = 'output/logs/020_norm/{library}-bbnorm.log'
    params:
        target = 60
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbduk_container
    shell:
        'bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq_norm} '
        'outt={output.fq_toss} '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target={params.target} '
        'min=5 '
        'peaks={output.peaks} '
        '2> {log.norm} '


# 01 trim and decontaminate reads
rule trim_decon:
    input:
        unpack(trim_resolver)
    output:
        fq = 'output/010_trim-decon/{library}.fq.gz',
        f_stats = 'output/010_trim-decon/{library}_filter-stats.txt',
        t_stats = 'output/010_trim-decon/{library}_trim-stats.txt'
    log:
        filter = 'output/logs/010_trim-decon/{library}-filter.log',  
        trim = 'output/logs/010_trim-decon/{library}-trim.log',
        reheader = 'output/logs/010_trim-decon/{library}-reheader.log'
    params:
        filter = '/phix174_ill.ref.fa.gz',
        trim = '/adapters.fa'
    threads:
        multiprocessing.cpu_count()
    singularity:
        bbduk_container
    shell:
        'reformat.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '2> {log.reheader} '
        '| '
        'sed -e \'s/^@\S*\s/@/g\' '          # remove SRA from fastq header
        '| '
        'reformat.sh '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'addcolon=t '
        'trimreaddescription=t '
        '2>> {log.reheader} '
        '| '
        'sed -e \'/^@\S\+:/s/$/N:0:NNNNNN/g\' '      # add second word to fastq header
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.f_stats} '
        '2> {log.filter} '
        '| '
        'bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'int=t '
        'out={output.fq} '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '

rule repair:
    input:
        r1 = 'output/fastq/{sample_name}/{sample_name}_1.fastq',
        r2 = 'output/fastq/{sample_name}/{sample_name}_2.fastq'
    output:
        r1 = 'output/fastq_repaired/{sample_name}_1.fastq.gz',
        r2 = 'output/fastq_repaired/{sample_name}_2.fastq.gz'
    threads:
        1
    resources:
        mem_gb = 50
    log:
        'output/logs/repair/{sample_name}.log'
    singularity:
        bbduk_container
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'zl=9 '
        'repair=t '
        '-Xmx{resources.mem_gb}g '
        '2> {log}'


rule dump_fastq:
    input:
        'output/SRAs/{sample_name}.sra'
    output:
        r1 = temp('output/fastq/{sample_name}/{sample_name}_1.fastq'),
        r2 = temp('output/fastq/{sample_name}/{sample_name}_2.fastq'),
        tmpdir = temp(directory('output/fastq/tmp_{sample_name}'))
    priority:
        1
    threads:
        multiprocessing.cpu_count()
    params:
        outdir = 'output/fastq/{sample_name}'
    log:
        'output/logs/dump_fastq/{sample_name}.log'
    singularity:
        sra_container
    shell:
        'fasterq-dump '
        '--outfile {wildcards.sample_name} '
        '--outdir {params.outdir} '
        '--temp {output.tmpdir} '
        '--threads {threads} '
        '--details '
        '--split-files '
        '--log-level 5 '
        '{input} '
        '&> {log} '

rule download_sra:
    output:
        temp('output/SRAs/{sample_name}.sra')
    params:
        url = lambda wildcards: name_to_url[wildcards.sample_name]
    threads:
        1
    log:
        'output/logs/download_sra/{sample_name}.log'
    shell:
        'wget '
        '-O {output} '
        '{params.url} '
        '&> {log}'
