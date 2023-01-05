# @author : johaGL

import os
from fun_tflex import *

SAMPLESBYMETADATA, speciesbymetadata = getsamplesNspe(config['metadata'])

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
del(SAMPLES)
del(RUN)
## -*- END TODO

assert len(set(speciesbymetadata)) == 1, "ERROR, you have more than 1 species in your samples metadata file"
species = speciesbymetadata[0]

# expected alignment files (output STAR)
EALF = ["Aligned.out.sam",  "Log.final.out",  "Log.out",  "Log.progress.out", "ReadsPerGene.out.tab",  "SJ.out.tab"]


rule all:
    input:
        expand(config['mappeddir']+"{sample_sp}/{ealf}", sample_sp=SAMPLESBYMETADATA, ealf=EALF) 

rule mapwithstar:
    input :
        t1 = config['trimmfqdir']+"{sample_sp}_trimm_"+PAIRS[0]+".fastq.gz",
        t2 = config['trimmfqdir']+"{sample_sp}_trimm_"+PAIRS[1]+".fastq.gz",
        gtf = config['gtf'] # config['gtf'][species]
    message :
        "  mapping  {wildcards.sample_sp}   ::"+species+"::"
    params :
        Idir = config['starindexdir'], #config['starindexdir'][species],
        prefix = config['mappeddir']+"{sample_sp}/",
        sthreads = config['star']['threads'],
        maxNbmap = config['star']['maxNbmap'],
        max_mismatches = config['star']['max_mismatches']
    output :
        sj = config['mappeddir']+"{sample_sp}/SJ.out.tab",
        al = config['mappeddir']+"{sample_sp}/Aligned.out.sam" ,
        lf = config['mappeddir']+"{sample_sp}/Log.final.out",
        lo = config['mappeddir']+"{sample_sp}/Log.out",
        lp = config['mappeddir']+"{sample_sp}/Log.progress.out",
        rp = config['mappeddir']+"{sample_sp}/ReadsPerGene.out.tab"
    shell :
          "STAR --runThreadN {params.sthreads} "
            "--quantMode GeneCounts "
             "--genomeDir {params.Idir} "
             "--sjdbGTFfile {input.gtf} "
             "--readFilesCommand gunzip -c "
             "--readFilesIn {input.t1},{input.t2} "
             "--outFilterType BySJout "
             "--outFilterMultimapNmax {params.maxNbmap} "
             "--alignSJoverhangMin 8 "
             "--alignSJDBoverhangMin 1 "
              "--outFilterMismatchNmax {params.max_mismatches} "
             "--outFilterMismatchNoverReadLmax 0.04 "
             "--alignIntronMin 20 "
             "--alignIntronMax 1000000 "
             "--alignMatesGapMax 1000000 "
             "--outFileNamePrefix {params.prefix}  "
