# tflex : a Snakemake pipeline for customized raw RNA-seq data pretreatment in cbib bb8 cluster

## Overview

Pipeline in active development, intended for RNA-seq analysis, in our laboratory.

Most of our datasets being paired-end, scripts have been tailored for them.

The objective is to allow flexibility in input file names, trying to keep them as much as possible
without the need to rename/relocate them. 

This pipeline must be kept into an intact directory, **DO NOT use it as a working directory**. To clone it : 
```
cd $HOME
git clone https://github.com/cbib/tflex
```
Alternatively, you can clone it into another folder, example let DIRSOFTWARE=/gpfs/home/myname/mysoftware/.

Usage : Go to your project `cd $MYPROJECT` and invoke a specific pipeline by typing `$DIRSOFTWARE/tflex/Snake_...smk`. The best is to do a dryrun before the real run, see below.

 
----------------------------------------------
# Developers README:

When two or more species are involved, conflicts arise. That is the reason why this snakemake is divided in parts. 
Each part is launched via bash scripts (making use of Slurm), see [cmd_examples/](cmd_examples/)

A yaml configuration file example is [configRat.yaml](https://github.com/johaGL/cocultureProj/blob/master/src1/configRat.yaml).

Software used in our pipelines: fastp, multiqc, STAR .
Versions are loaded with `module load` in .sh final scripts, by user. For example, STAR versions 2.7.1a and 2.7.9a have already been tested with success.

**TODO**: I will provide an example directory inhere, with a simplistic suggested usage. 

## Pipelines, their rules and their options:
1. `Snake_doindex.smk` : needs a .yaml file exclusively for Index build-up. That .yaml must contain URLs for fasta, gtf and gff3 files.

2. `Snake_fastp.smk` : options that pass from .yaml to fastp options:
  * threads : `-w` which means workers
  * minqual : `--qualified_quality_phred` and `--cut_mean_quality`
  
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

| strucfile | explanation | 
| -------  | ---------- | 
| 1 | if {sampleid}\_{pair name}.fastq.gz examples: `D123546_1.fastq.gz` , `T9_876-sansglu_R1.fastq.gz` |
| 2 | if {sampleid}\_{pair name}\_{run}.fastq.gz example: `T9_876-sansglu_R1_001.fastq.gz` , **where 001 is {run}** |

  In this way, {pair name} is detected automatically (can be 'R2' or '2' or equivalent).
   Only tested for  {run} being unique, up to now. If two or more {run} values, a new option must be created and configured, or simply use separated run folders. 
 
For full instructions about setting a config file see *Users README* section, "Micro-tutorial" subsection.

### Running a dryrun
User is encouraged to perform a dryrun, as it allows to detect errors. 




-----------------------------------------------

## Users README:


## * Part I * : QC, Mapping, Raw Counts Matrix

## I.1. Prepare your work

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
But improved (2023) `.sh` files are available at the end of this readme.

Following steps will help you to prepare and obtain a folder structure as shown above.
  * A. Create and save your metadata csv file(s)
  * B. Create the yaml configuration file(s):
      - for index building if not done.
      - for all the other (qc, mapping, counts) a unique yaml.    
      Example : configIndexHuman.yaml, configHuman.yaml. 
      Do separate files for each species (mouse, zebrafish, rat, etc).  
  * C. Create a 'dryrun' bash script.
  * D. Create your final bash script(s).

### A. create and save your metadata csv files
Files must be comma delimited, with a header. 
Here a file samples_mouse.csv:
```
sample,species,condition
DT44441,mouse,control
DT44442,mouse,control
DT44443,mouse,starved
DT44444,mouse,starved
```
Here for samples_human.csv:
```
sample,species,condition
DT55551,human,control
DT55552,human,control
DT55553,human,starved
DT55554,human,starved
```

### B. the configuration yaml file(s)
* a) for the Index (if not already available), example : [config_doindex_rat.yaml](https://github.com/johaGL/cocultureProj/blob/master/src1/config_doindex_rat.yaml). You must set URLs to ftp ensembl location of desired genomes (fasta) and transcripts coordinates (gff3, gtf). You can also use ENCODE as source.
* b) for QC, mapping and counts: a single .yaml must be created, one by species. Follow micro-tutorial below.
An example for b) [configRat.yaml](https://github.com/johaGL/cocultureProj/blob/master/src1/configRat.yaml)

#### Micro-tutorial to set your yaml file for QC+mapping+count (estimated reading time 30 seconds):

Use this template for your config_xxx.yaml
```
strucfile :  
fastqdir : 
trimmfqdir : 
qc_reportsdir : 
mappeddir : 
fastp :
    threads : 
    minqual : 
star :
    threads : 
    maxNbmap : 
    max_mismatches : 999

# part species specific : 
metadata :  
starindexdir :  
gtf : 
```
In your config_xxx.yaml, set your options : 

  * strucfile

 | strucfile | explanation | 
| -------  | ---------- | 
| 1 | if {sampleid}\_{pair name}.fastq.gz examples: `D123546_1.fastq.gz` , `T9_876-sansglu_R1.fastq.gz` |
 | 2 | if {sampleid}\_{pair name}\_{run}.fastq.gz example: `T9_876-sansglu_R1_001.fastq.gz` , **where 001 is {run}** |
 
  * fastp 
     - threads : `-w` which means workers, for fastp
     - minqual : `--qualified_quality_phred` and `--cut_mean_quality`, for fastp
  * star
     - threads : `--runThreadN`, for STAR
     - maxNbmap : `--outFilterMultimapNmax`, for STAR
     - max_mismatches : `--outFilterMismatchNmax`, for STAR
  * fastqdir, trimmfqdir, qc_reportsdir, mappeddir, indexdir must be directory names.
  * metadata and gtf must be file names, in absolute path.

 When max_mismatches set to 999, it  deactivates --outFilterMismatchNmax, so by default 4% mismatches are accepted (because Snake_starmp.smk fixed  --outFilterMismatchNoverReadLmax 0.04) .
 
Remember:  
1. Respect indentations
2. Make sure that every directory name ends with "slash" `/`
3. Integers are not quoted
4. Strings ARE quoted in .yaml file
5. All files and directories must be written in **absolute** paths
6. Make sure output folders either absent or empty
7. If you feel lost see example  [configRat.yaml](https://github.com/johaGL/cocultureProj/blob/master/src1/configRat.yaml)


### C.  create your 'dryrun' script
This is **STRONGLY** recomended, as it will not run for real 
but it will show the steps that will be executed and detect any problems. 
As no full calculations will be executed, it does not need sbatch.
Example: [dryrunSteps1_2.sh](cmd_examples/dryrunSteps1_2.sh) 

### D. write definitive script(s)
Example : [step1.sh](cmd_examples/step1.sh)

Adapt SLURM options to your project. Use same commands as written in your dryrun, make sure your dryrun is clean and that output folders are absent or empty (our pipelines do not overwrite!).


## I.2. Launch definitive script
First make sure all steps above have been sucessful, including the dryrun ! 
If no **dryrun** errors, congrats! Make sure the **output** folders are empty, because the pipelines will not overwrite, and just go for a cup of coffee after launching `sbatch step1.sh` .  


## * Part II * : Downstream analysis (in development)

create the conda environment : 
```
conda env create -f r_env_tflex.yml
```

For all the downstream analysis, see [Rscripts/](https://github.com/cbib/tflex/blob/master/Rscripts/)).


1. In your project folder, create scr2/ folder and put inside:
 - the .yml file(s) such the one shown at the end of this section "Example yml downstream". Each performs one comparison (or seeks for the existing dds object that has the comparison inside, by using the optional parameter `filehasthiscontrast` )
 
2. In cbib cluster activate the `r_env_tflex` environment. 

3. Launch by comparison: 
```
Rscript $TFLEX_SCRIPT $MY_YML $LOCATION_TFLEX 
```

Use the scripts in Rscripts in sequential order. 


"Example yml downstream"

Here two different yml, each for one comparison.

- the dds object does not exists, I launch the `comparison_MvsV.yml` : 
```
mywdir : "~/bulk202303_TD/myelin_202303/" 
metadatafile : "raw_counts/metadata_myelin.csv" 
countsfile : "raw_counts/rawcounts_myelin.tsv"

conditionLEVELS : [ "Vehicle", "Pellets" , "Myelin" ] # "control" here first
requiredcontrast : ['condition','Myelin', 'Vehicle' ] # control goes last
outname :  "mpv"
shortname : "human"
SPECIESensmbldsetname : "hsapiens_gene_ensembl"
equivalentid : "ensembl_gene_id_version"
# note: if ensembl ids have a dot '.' equivalentid : "ensembl_.._version"

samplestodrop : 
  - "Vehicle_24h-3"

nb_bt : 100 #min nb of genes req for a biotype, for biotype to be displayed

### vars  in 2_ddstotabplots.R
trustedcolumn : "symbol_unique" 
# for tables "top" (tsv files):
absLFCcutoff : 0
sigcutoff : 0.05

# cutoffs for volcano labels
aLFClab : 0.5
pjlab : 0.05
labelsinvolcano : True

# PCA related 
maxncontrib : 5
#cutoffsvarspvals : [0.05, 0.01, 0.005,0.001, 0.0005, 0.0001]
dobiplots_bypval_bycond : False
```

- the dds object already exists (contrast is different but contained in the previously generated dds) `comparison_PvsV.yml`:
```
mywdir : "~/bulk202303_TD/myelin_202303/" 
metadatafile : "raw_counts/metadata_myelin.csv" 
countsfile : "raw_counts/rawcounts_myelin.tsv"

conditionLEVELS : [ "Vehicle", "Pellets" , "Myelin" ] # "control" here first
requiredcontrast : ['condition','Pellets', 'Vehicle' ] # control goes last
objecthasthiscontrast : results_mpv/rds/ddsObj_Myelin_vs_Vehicle_mpv.rds
outname :  "mpv"
shortname : "human"
SPECIESensmbldsetname : "hsapiens_gene_ensembl"
equivalentid : "ensembl_gene_id_version"
# note: if ensembl ids have a dot '.' equivalentid : "ensembl_.._version"

nb_bt : 100 #min nb of genes req for a biotype, for biotype to be displayed

### vars  in 2_ddstotabplots.R
trustedcolumn : "symbol_unique" 
# for tables "top" (tsv files):
absLFCcutoff : 0
sigcutoff : 0.05

```

--------------------------------------
#### Improved .sh files (with Slurm)  for steps QC, trimming, mapping

- [step1](cmd_examples/step1.sh)
- [step2](cmd_examples/step2.sh)
- [step3](cmd_examples/step3.sh)


----------------------------------------
Acknowledgments:

to the creators of snakemake.
to the authors of cited software.

--------------------------------------------
####  Authors
Joha GL , B Dartigues.

Many thanks to A Barré for help with Slurm.


