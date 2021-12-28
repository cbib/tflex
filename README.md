# Snakemake pipeline for customized raw RNA-seq data pretreatment in cbib cluster

## Overview

Pipeline in active development, intended for RNA-seq fastq files quality control, mapping,
possible to extend to other applications.

Most of our datasets being paired-end, scripts have been tailored for them.

The objective is to allow flexibility in input file names, trying to keep them as much as possible
without the need to rename them. 

This pipeline must be kept in an intact directory, do not use it as a working directory.

Working directory must be a RNA-seq analysis project located
completely outside **tflex**. Working directory can contain 
your data and/or results, in any case, all paths must be defined in a configuration `.yaml` file.
An example of such "project" (is simply a directory) can be seen `here`. 
In other words, lets say your working directory is '$HOME/aflyproject/' you invoke this tflex for example in this way:
```
$ cd $HOME/
$ git clone .../tflex.git
$ ls $HOME/aflyproject
$ >     config_doindex.yaml  config_qcmap.yaml  fastqfiles/
$ cd $HOME/aflyproject
$ snakemake -s ../tflex/Snake_doindex --configfile config_doindex.yaml --cores 1' 
```

## tflex

You are here in tflex, it does not contain any data, user files nor configuration files
As explained above, DO NOT ADD ANY FILES, unless you know what you are doing.



###  




#### ---- Authors
Joha GL 
... B Dartigues ...  A Barré ...
