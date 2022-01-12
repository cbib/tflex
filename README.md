# Snakemake pipeline for customized raw RNA-seq data pretreatment in cbib bb8 cluster

## Overview

Pipeline in active development, intended for RNA-seq fastq files quality control, mapping,
possible to extend to other applications.

Most of our datasets being paired-end, scripts have been tailored for them.

The objective is to allow flexibility in input file names, trying to keep them as much as possible
without the need to rename them. 

This pipeline must be kept into an intact directory, **DO NOT use it as a working directory**. To clone it : 
```
cd $HOME
git clone https://github.com/johaGL/tflex
```
Alternatively, you can clone it into another folder, example let DIRSOFTWARE=/gpfs/home/myname/mysoftware/, then go to your project (cd \$MYPROJECT) and to invoke a specific pipeline you must type `$DIRSOFTWARE/tflex/Snake_...smk`.
 
**Note for users** : tflex must only be composed of **tools**, but NOT fastq files or any other files proper to your RNA-seq project. For specific needs, fork it and feel free to modify/contribute/suggest.  

# Developers README:

Pipelines are designed to  be used into a bash script like [this one](https://github.com/johaGL/cocultureProj/blob/e0418b945e01ab9cf0165f41b28d64c35091be34/src1/step1.sh).

A yaml configuration file example is [configRat.yaml](https://github.com/johaGL/cocultureProj/blob/master/src1/configRat.yaml).

Software used in our pipelines: fastp, multiqc, STAR .
Versions are loaded with `module load` in .sh final scripts, by user. For example, STAR versions 2.7.1a and 2.7.9a have already been tested with success.

**TODO**: I will provide an example directory inhere, with a simplistic suggested usage. 

## Pipelines, their rules and their options:
1. `Snake_doindex.smk`
2. `Snake_fastp.smk` : options that pass from .yaml to fastp options:
  * thr_p : `-w` which means workers
  * minqual : `--qualified_quality_phred` and `----cut_mean_quality`
  
3. `Snake_starmp.smk` :  options that pass from .yaml to STAR options:
  * threads : `--runThreadN`
  * maxNbmap : `--outFilterMultimapNmax`
  * max_mismatches : `--outFilterMismatchNmax`

4. `countsfromstarout.R` This is Rscript **TODO** integrate to snakemake 


## A  config.yaml template
```
strucfile : 1 
fastqdir : ""
trimmfqdir : ""
qc_reportsdir : ""
mappeddir : ""
fastp :
    threads : 12
    minqual : 30
star :
    threads : 12
    maxNbmap : 1
    max_mismatches : 999

# part species specific : 
metadata :  ""
starindexdir:  ""
gtf : ""
```
 

### Comments about this yaml template:
* strucfile : is an integer, to set the option about the files names "structure", up to now I have defined two possible options:

   1 --> if {sampleid}\_{pair name}.fastq.gz
	  examples: `D123546_1.fastq.gz` , `T9_876-sansglu_R1.fastq.gz`
   2 --> if {sampleid}\_{pair name}_{run}.fastq.gz
	  example: `T9_876-sansglu_R1_001.fastq.gz` , **where 001 is {run}**

  In this way, {pair name} is detected automatically (can be 'R2' or '2' or equivalent).
   Only tested {run} unique, up to now. If two or more "run" values, a new option must be created and configured, or simply use separated run folders. 
   
* User must guarantee every directory name ends with `/`

 * In "star" key, max_mismatches set to 999  deactivates --outFilterMismatchNmax, so by default 4% mismatches are accepted


# Users README:

## Prepare your work

* User working directory must be a RNA-seq analysis project located completely outside tflex.
* Working directory can contain your data and/or results, there is not a mandatory folder structure because in any case, all paths must be defined in configuration `.yaml` file(s).

Lets say your working directory is MYPROJECT=/gpfs/home/myname/myproject with an organization as follows : 

```
$ tree -L 2 $MYPROJECT

├── metadata
│   ├── samplesbothspecies.csv
│   ├── samplesspecies1.csv
│   └── samplesspecies2.csv
├── README.md
├── referencegenomes
│   └── Rattus_norvegicus.mRatBN7.2
└── src1
    ├── cmd_doindex_rat.sh
    ├── config_doindex_rat.yaml
    ├── configHuman.yaml
    ├── configRat.yaml
    ├── dryrunSteps1_2.sh
    ├── slurm_output
    ├── step1.sh
    └── step2.sh
```

An example of such "project" can be seen [here](https://github.com/johaGL/cocultureProj/). 
Follow this steps to prepare and obtain a folder structure as shown above.
  * A. Create and save your metadata
  * B. Create the yaml configuration files:
      - for index building if not done.
      - for all the other (qc, mapping, counts) a unique yaml.    
      Example : configIndexHuman.yaml, configHuman.yaml. 
      Do separate files for each species (mouse, zebrafish, rat, etc).  
  * C. Create a 'dryrun' bash script.
  * D. Create your final bash script(s).

### A. create and save your metadata csv files
Files must be comma delimited, with a header. 
Here a file for one species:
```
sample,species,condition
DT44441,mouse,control
DT44442,mouse,control
DT44443,mouse,starved
DT44444,mouse,starved
```
Here for the other one:
```
sample,species,condition
DT55551,human,control
DT55552,human,control
DT55553,human,starved
DT55554,human,starved
```

### B. the configuration yaml file(s)
* a) for the Index (if not already available)
* b) for QC, mapping and counts: a single .yaml must be created, one by species.
An example for b) [configRat.yaml](https://github.com/johaGL/cocultureProj/blob/master/src1/configRat.yaml)

IMPORTANT: make sure that every directory name ends with "slash" `/`

### C.  create your 'dryrun' script
This is **STRONGLY** suggested, as it will not run for real 
but it will show the steps that will be executed and detect any problems. 
As no real calculations will be executed, it does not need sbatch.
Example: [dryrunSteps1_2.sh](https://github.com/johaGL/cocultureProj/blob/master/src1/dryrunSteps1_2.sh) 

### D. write definitive script(s)
Example : [step1.sh](https://github.com/johaGL/cocultureProj/blob/e0418b945e01ab9cf0165f41b28d64c35091be34/src1/step1.sh)

Adapt SLURM options to your project.



## Launch definitive script
First make sure all steps above have been sucessful, including the dryrun ! 
If no errors, congrats! Make sure the **output** folders are empty, because the pipelines will not overwrite, and just go for a cup of coffee after launching `sbatch step1.sh` .  

--------------------------------------------
####  Authors
Joha GL , B Dartigues, A Barré.
