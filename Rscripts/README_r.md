# Downstream analysis scripts : tflex/Rscripts/ folder (routine phase)


Keeping in line with our example:

* `cocultureProj` is the place where:
    + Count matrices (`datarc/`), and metadata must be present (code is currently designed to accept relative paths for these folders, but in the future, think about changing scripts to make them accept absolute paths).
    + scripts, gmt files and biomart output files are saved (in their respective sub-folders: `scr(...)/`, `gmt_saved/`, `biomart_bm/`). File format .rds is preferred because it offers a level of compression. The gmt files are not sync (.gitignore) owing to their large size. 
    + all `results_` pre-fixed folders and their nested folders will be saved. 
* by default labels in volcano are plotted. If this behaviour is not desired, set `labelsinvolcano` to `False`.
* All automatic Rscripts accept same yaml file.


```
cocultureProj
├── biomart_db
│   ├── bm_human.rds
│   ...
│   └── orthoinfo_rat.rds
├── datarc
|   ...
│   ├── raw_countsHuman-new.tsv
│   └── raw_countsRat-new.tsv
├── gmt_saved
│   ├── gmtfiltered_CP:KEGG-v2-human.tsv
|   ...
│   └── info-filter_GO:BP-v2-rat.txt
├── metadata
|   ...
│   ├── samplesRat.csv
│   ├── samplesRat-new.csv
│   └── table_antonio.csv
├── README.md
├── results_Neurons_rat
│   ├── plots
│   ├── rds
│   └── tables
├── results_P3_human-CD3s
│   ├── plots
│   ├── rds
│   └── tables
...
├── scr2     <=============== * * * you are here * * *
│   ├── 0_verify_options.R
│   ├── 1_dataprep.R
│   ├── 2_a_DE.R
|   ...
│   ├── func.R
│   ├── misc_scr2
│   ├── NAU0_rename.R
...
│   ├── README_scr2.md
│   ├── set1_launch.sh
│   ├── set1_Neu_insert-vs-control.yaml
│   ├── set1_P3-CD_INSERT-vs-controlIns.yaml
...
│   └── yamltemplate
...
```

## Automated Rscripts :
* `0_verify_options.R` : prints the options (required) that you defined in your configuration file.
* `1_dataprep.R` : connects ensembl id to biotypes via biomaRt (biomart_bd folder is locally kept). Saves plots (horizontal bars, multivariate initial on rlog), 
* `2_a_DE.R` runs DE and saves DESeq2 object (rds).
* `2_b_tabplots.R` : uses DESeq2 saved object to perform a volcano plot, saves tables.
* `3_a_PCAonNormalized.R` : PCA on cpm (counts per million) matrix
* `3_b_plotPCA.R`
*  all other starting by a number n `n_`.

These automatic Rscripts constitute a sequence of steps, using a YAML configuration file. All automatic Rscripts accept same yaml file.

Configuration **.yaml** file is mandatory. Here an example:
```
mywdir : "~/cocultureProj/" 
metadatafile : "metadata/samplesRat.csv"
countsfile : "datarc/raw_countsRat.tsv"
conditionLEVELS : ['control' , 'insert'] # control is first, to order factors
requiredcontrast : ['condition', 'insert', 'control'] # control goes last
outname :  "Neurons_rat"
shortname : "rat"
SPECIESensmbldsetname : "rnorvegicus_gene_ensembl"
equivalentid : "ensembl_gene_id"
# note: if ensembl ids have a dot '.' equivalentid : "ensembl_.._version"

nb_bt : 100 #min nb of genes req for a biotype, for biotype to be displayed

### vars  in 2_ddstotabplots.R
trustedcolumn : "symbol_unique" 
absLFCcutoff : 1.1
sigcutoff : 0.05
aLFClab : 1.5
pjlab : 4e-04
labelsinvolcano : True
```


Example:
```
$ Rscript --vanilla tflex/Rscripts/1_*.R  $YAMLFILE
```
Once 1\_ finishes, if everything went ok you can launch 2, 3 (sequentially):
```
$ Rscript --vanilla tflex/Rscripts/2_*.R  $YAMLFILE
```






