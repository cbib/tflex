

import os
#import pandas as pd
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


species1 = config['species'][0]
species2 = None
if len(config['species']) == 2: # this version supports max 2 species
    species2 = config['species'][1]

#xixo = pd.read_csv(config['metadata'])
#print(xixo.shape)
#print("======")
dicosp = d2o(config['metadata'], config['species'])

SPECIES1 = dicosp[species1] # list of samples this species
SPECIES2 = dicosp[species2] # list of samples this other species
EOSAM = ["Aligned.out.sam",  "Log.final.out",  "Log.out",  "Log.progress.out", "ReadsPerGene.out.tab",  "SJ.out.tab"]
def prepmapinfo(x,y, config):
   """x : dictionary,  y : species, config : configuration dictionary"""
   return "Preparing for Mapping, the number of samples is: "+str(len(x[y]))+\
            "; species is: "+y+";  gtf: "+config['gtf'][y]+'\n'+\
          "STARindex: "+config["starindexdir"][y]

sp1message = prepmapinfo(dicosp, species1, config)
sp2message = prepmapinfo(dicosp, species2, config)


##########  RULES

rule all:
    input:
        # expand(config['trimmfqdir']+"{sample}_trimm_{pair}.fastq.gz", sample=SAMPLES, pair=PAIRS),
        expand(config['qc_reportsdir']+"{sample}/fastp.{exten}", sample=SAMPLES, exten=['html','json']),
        config['qc_reportsdir']+"multiqc_fastp/multiqc_report.html",
        expand(config['trimmfqdir']+"{sa_spe1}_trimm_{pair}.fastq.gz", sa_spe1=SPECIES1, pair=PAIRS), 
        expand(config['prep4mapdir']+"infosspecies-{sp}.txt", sp=config['species']),
        
        expand(config['mappeddir']+"{sa_spe1}/{eosam}", sa_spe1=SPECIES1, eosam=["SJ.out.tab"] ),  # !!!!
        config['prep4mapdir']+"mappping-"+species1+"-ended.txt",
        expand(config['trimmfqdir']+"{sa_spe2}_trimm_{pair}.fastq.gz", sa_spe2=SPECIES2, pair=PAIRS)
# TODO: ,expand(config['mappeddir']+"{sa_spe2}/Aligned.out.sam", sa_spe1=SPECIES2)
        
       
# https://stackoverflow.com/questions/57028398/how-can-i-run-a-subset-of-my-snakemake-rules-several-times-with-wildcards

rule fastp:
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


rule checkpoint_species1:
    input : 
       qcreport = rules.multiqc.output.qcreport
    message : sp1message
    output:
        finf = config['prep4mapdir']+"infosspecies-"+species1+".txt"
    params:
        sms = sp1message+'\nSamples:\n'+'\n'.join(dicosp[species1])
    shell: "printf '{params.sms}\n' > {output.finf}"

rule mapping_species1:
    input  :   
        t1 = config['trimmfqdir']+"{sa_spe1}_trimm_"+PAIRS[0]+".fastq.gz",
        t2 = config['trimmfqdir']+"{sa_spe1}_trimm_"+PAIRS[1]+".fastq.gz",
        gtf = config['gtf'][species1] 
    message : 
        "  Mapping  {wildcards.sa_spe1}   ::"+species1+"::" 
    output : 
        sjo = config['mappeddir']+"{sa_spe1}/SJ.out.tab",
        #ad = config['prep4mapdir']+"mappping-"+species1+"-ended.txt"
        fiu = config['prep4mapdir']+"mappping-rat-ended.txt" 
    params :
        Idir = config['starindexdir'][species1],
        prefix = config['mappeddir']+"{sa_spe1}/",
        sthreads = config['star']['threads'],
        max_mismatches = config['star']['max_mismatches']

    
    shell : 
          " touch {output.sjo}"
