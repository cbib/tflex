# @author : johaGL

import os
from fun_tflex import *

# 'strucfile' option must be specified in configfile:
#     note :  a {pair name} example is 'R2' or '2' or equivalent
#   1 -- > if {sampleid}_{pair name}.fastq.gz
#	  examples: D123546_1.fastq.gz , T9_876-sansglu_R1.fastq.gz
#   2 -- > if {sampleid}_{pair name}_{run}.fastq.gz
#	  example: T9_876-sansglu_R1_001.fastq.gz , so {run} is '001'

# REMARK: rules have been tested only if 'run' is the same across  all files, or absent
   

# PULL wildcards function
def pullwildcardsbyopt(strucfile, rawdir):
    """strucfile is fastq file pattern, an integer (choice)"""
    if strucfile == 1:
        FASTQNAMES, PAIRNAMES, = glob_wildcards(os.path.join(rawdir,'{saname}_{pair}.fastq.gz'))
        RL = [ '' ] # empty string in run list means absent {run} in filename
    elif strucfile == 2:
        FASTQNAMES, PAIRNAMES, RL, = glob_wildcards(os.path.join(rawdir, '{saname}_{pair}_{run}.fastq.gz'))
    return makelistsunique(FASTQNAMES, PAIRNAMES, RL)

########## global variables (from configuration file)

strucfile = config['strucfile']
rawdir = config['fastqdir']

SAMPLES, PAIRS, RUN, = pullwildcardsbyopt(strucfile, rawdir)
assert len(set(RUN)) == 1, "CHECK CODE ! {run} suffix is not a constant, never tested before"
assert len(sorted(set(PAIRS))) == 2, "Error!! PAIRS list must be length 2, and ordered"
RUN = listtostr_single(strucfile, RUN) # string, function has assert part


rule all:
    input:
        expand(config['trimmfqdir']+"{sample}_trimm_{pair}.fastq.gz", sample=SAMPLES, pair=PAIRS),
        expand(config['qc_reportsdir']+"{sample}/fastp.{exten}", sample=SAMPLES, exten=['html','json']),
        config['qc_reportsdir']+"multiqc_fastp/multiqc_report.html"

rule fastp:
   """improvement: added cut_front and cut_tail"""
    input :
        reads1 = rawdir+"{sample}_"+PAIRS[0]+RUN+".fastq.gz",
        reads2 = rawdir+"{sample}_"+PAIRS[1]+RUN+".fastq.gz"
    output:
        trimm1 = config['trimmfqdir']+"{sample}_trimm_"+PAIRS[0]+".fastq.gz",
        trimm2 = config['trimmfqdir']+"{sample}_trimm_"+PAIRS[1]+".fastq.gz",
        ohtml = config['qc_reportsdir']+"{sample}/fastp.html",
        ojson = config['qc_reportsdir']+"{sample}/fastp.json"
    message :
        "    **    processing {wildcards.sample}    **    "
    params:
        thr_p = config['fastp']['threads'],
        minqual = config['fastp']['minqual']
    shell:
        "fastp --detect_adapter_for_pe -w {params.thr_p} \
              --qualified_quality_phred {params.minqual} \
              --cut_front --cut_tail --cut_mean_quality {params.minqual} \
              -i {input.reads1} -I {input.reads2} \
              -o {output.trimm1} -O {output.trimm2} \
              --html {output.ohtml} --json {output.ojson}"

rule multiqc:
    input : 
        trimmedfiles = expand(config['trimmfqdir']+"{sample}_trimm_{pair}.fastq.gz", sample=SAMPLES, pair=PAIRS)
    output:
        qcreport = config['qc_reportsdir']+"multiqc_fastp/multiqc_report.html"
    params:
        indir = config['qc_reportsdir'],
        odir = config['qc_reportsdir']+"multiqc_fastp/"
    shell:
        "multiqc -f --outdir {params.odir} {params.indir}"
