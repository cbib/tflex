# Bunch of functions for handling counts and associated movements
# johaGL
# set the class
setClass(Class="rnaexp", slots=list(rawq="array", 
                                    metadata="data.frame",
                                    row_data="data.frame",
                                    cpm="array"))

o_getcols <- function(myobj, ...){
  cat("\nusage : o_getcols(obj, vectorofcols) # yields cols in vectorofcols ",
      "\n    o_getcols(obj,'-',vectorofcols) # yields cols NOTin vectorofcols \n")
  argus <- list(...)
  if(length(argus) ==  1 ){ 
    newmeta <- myobj@metadata %>% filter(sample %in% argus[[1]])
  }else if (length(argus) == 2 ){ # a minus was entered
    if(argus[[1]][1] =="-"){
      newmeta <- myobj@metadata %>% filter(! sample %in% argus[[2]]) #!
    }else{"Error! if entering three arguments, the 2nd must be '-', see usage" }
  }else{ print("ERROR ! o_getcolumns  accepts max 3 arguments ") } 
  newmat <- myobj@rawq[,newmeta$sample]
  return(new("rnaexp", rawq=as.array(newmat, dim=dim(newmat)), 
             metadata=newmeta,
             row_data=myobj@row_data))
}

o_getrows <- function(myobj, criterion, myvector){
  # for the moment assuming rawq has symbol_unique as rownames
  cat("\nusage : o_getrows(obj,'symbol_unique', c('BRCA2', ...))", 
      "\n    o_getrows(obj, 'ensembl_gene_id_version', c('ENSG..9.23',...))",
      "\n    o_getrows(obj, 'biotype', c('protein_coding') \n")
  
  if(criterion !=  "biotype" ){ 
    newrd <- myobj@row_data[match( myvector,myobj@row_data[[criterion]]),]
  }else if (criterion == "biotype"){ 
    newrd  <- myobj@row_data %>% filter(gene_biotype %in% myvector) 
  }
  newmat <- myobj@rawq[newrd$symbol_unique,]
  rownames(myobj@row_data) <- NULL
  return(new("rnaexp", rawq=as.array(newmat, dim=dim(newmat)), 
             metadata=myobj@metadata,
             row_data=newrd))
}

corner <- function(amatrix){
  print(amatrix[1:5,1:5])
}

savetsv <- function(elem, namefi){
  write.table(elem, namefi, sep='\t', col.names=TRUE, row.names=FALSE )
}
opentsv <- function(namefi){
  print("set check.names=F in read.table to avoid replacing '-' by '.' ")
  return(read.table(namefi, sep='\t', header=TRUE, check.names = FALSE))
}

lastnsubstr <- function(mystr, n){
  # lastnsubstr(mystr,n) yields the last n chars of mystring, uses base::substr
  newx <- substr(mystr, nchar(mystr) - n + 1, nchar(mystr)) 
  return(newx)
}

getoutputdirs <- function(outer=o_dir, inner=c("rds/","tables/", "plots/")){
  # create/get full paths vector in this order: rds/ tables/ plots/
  dirs_ <- sapply(inner, function(x) paste0(o_dir, x))
  for (d in dirs_){
    if (!dir.exists(d)){
      dir.create(d, recursive = TRUE)
    }
  }
  return(dirs_)
}

print_warnings_rawcdfANDmetadata <- function(rawc, metadata){
  # to explore conveniently rows and cols in *counts dataframe* and metadata
  dim(rawc)
  head(rawc)
  print("* Exploring ensembl_id column in raw counts (uniqueness, etc)*")
  if (all(rawc$ensembl_id==unique(rawc$ensembl_id))){
    print("ok ensembl_id are unique in raw counts matrix")
  }else{print("WARNING! : ensembl_id are NOT unique")}
  print(colnames(rawc))
  ids <- tibble(id = rawc$ensembl_id)
  withpoints <- ids %>% filter(str_detect(id,"\\."))
  if (nrow(withpoints) == 0){
    print("NOTdetected dot '.' in id, set equivalentid : 'ensembl_gene_id' in your .yaml")
  }else{print("detected dot '.' in id, set equivalentid : 'ensembl_gene_id_version' in your .yaml")}
  if (!all(metadata$sample %in% colnames(rawc))){
    print("ERROR, sample names in metadata and rawcounts do not match ")}
}

countsdftomatrix <- function(countsdf){
  # df to true matrix, NA replaced by zero, all-zero rows excluded.
  print("countsdftomatrix RETURNS: matrix, and drops 'zeroes-all rows'")
  print(paste("function countdftomatrix assumes:",
              "first column contains rownames and the rest is numeric"))
  raw.mat <- countsdf
  rownames(raw.mat) <- raw.mat[,1]
  raw.mat <- raw.mat[, -1] 
  raw.mat[is.na(raw.mat)] <- 0 # replace NA values
  tokeep <- apply(raw.mat, 1, function(x) sum(x) > 0)
  raw.mat <- raw.mat[tokeep, ] # reject zeroes rows
  print(dim(raw.mat))
  outmat <- array(unlist(raw.mat), dim=dim(raw.mat),
               dimnames=list(rownames(raw.mat),colnames(raw.mat)))
  return(outmat)
}

getfullRowsSymbols <- function(rawmat, equi_id, bm_sp){
  cat("\ngetfullRowsSymbols : returns dataframe same order as matrix rownames")
  cat("\nadds column 'symbol_unique', but if absent symbol, value is ", equi_id,'\n')
  # note: equi_id can be ensembl_gene_id or ens_..id_version, depends.
  mybmdf <- data.frame(x = rownames(rawmat))
  mybmdf[[equi_id]] <- mybmdf$x
  print(colnames(bm_sp))
  print(head(mybmdf))
  mybmdf <- left_join(mybmdf, bm_sp, keep = FALSE) %>% distinct(.keep_all = T)
  vectorofsymbols <- c() ; for (i in 1:length(mybmdf[[equi_id]])){
    if (mybmdf[i,]$external_gene_name == '' | 
        is.na(mybmdf[i,]$external_gene_name)){
      vectorofsymbols <- c(vectorofsymbols, mybmdf[i,][[equi_id]])
    }else{ vectorofsymbols <- c(vectorofsymbols, mybmdf[i,]$external_gene_name)}
  }
  mybmdf$symbol_unique <- make.unique(vectorofsymbols)
  print("finished: column x ($x) has the common elements with the matrix")
  return(mybmdf)
}

furtherRowsSymbols <- function(rowdatadf, equi_id, bm_sp){
  print(paste("function getfullRowsSymbols() trusts in ",equi_id, "because unique",
              "but can miss things : ENSG00000179314.15==NA (ENSG00000179314.16==WSCD1) "))
  print("searching further with ensembl_gene_id only")
  nono <- rowdatadf %>% filter(is.na(external_gene_name) | external_gene_name == '')
  savedesperate <- c(); for (i in 1:dim(nono)[1]){
    origname = nono[i,equi_id]
    GENEIDONLY = unname(unlist(str_split(origname, '\\.'))[1])
    goodline <- bm_sp[bm_sp$ensembl_gene_id == GENEIDONLY,]
    if (dim(goodline)[1] > 0 ){
      colstofill <- c("ensembl_gene_id", "entrezgene_id",
                      "external_gene_name", "description", "gene_biotype")
      for (col in colstofill){
        rowdatadf[rowdatadf[[equi_id]] == origname, col ] <- goodline[[col]]
      }
      savedesperate <- c(savedesperate, origname)
    }
  } # end for
  # column symbol_unique : fill gaps and make unique
  rowdatadf <- rowdatadf %>% mutate(symbol_unique= case_when(
    (external_gene_name == "" | is.na(external_gene_name)) ~ symbol_unique,
    TRUE ~ external_gene_name
  ))
  rowdatadf$symbol_unique <- make.unique(rowdatadf$symbol_unique)
  return(list(rowdatadf, savedesperate))
}

do_biotype_plot <- function(gene_biotype, nb_bt, mytitle){
  # rawmat, bm_sp, mytitle, nb_bt, equi_id
  
  # from count matrix,bm outputs ggplot2 object. nb_bt :nb genes in biotype
  # mybmdf <- data.frame(x = rownames(rawmat))
  # mybmdf[[equi_id]] <- mybmdf$x
  # print(colnames(bm_sp))
  # print(dim(mybmdf))
  # print(dim(rawmat))
  # mybmdf <- left_join(mybmdf, bm_sp, keep = FALSE) %>% distinct(.keep_all = T)
  
  biotyfreq <- table(gene_biotype)
  sortedfq <- sort(biotyfreq[biotyfreq >= nb_bt])
  bty <- data.frame("Biotype"=factor(names(sortedfq), levels=names(sortedfq)),
                    "Count"=as.vector(sortedfq)) 
  
  if(nrow(bty)>=3){
    somecols <- rev(brewer.pal(nrow(bty),"Paired"))
    gp <- ggplot(bty)  + 
      geom_bar(aes(x=Biotype, y=Count, 
                   fill=somecols), stat="identity") +
      coord_flip() + theme_light() + theme(legend.position = "none") +
      labs(title= mytitle)
  }else{
    gp <- ggplot(bty)  + 
      geom_bar(aes(x=Biotype, y=Count, fill="darkcyan"), stat="identity") +
      coord_flip() + theme_light() + theme(legend.position = "none") +
      labs(title= mytitle)
  }
  return(gp)
}

symbomore4DE <- function(res.here, myobj, trustedcolumn, equi_id, 
                         exclude_types, tabl_dir, outname){
  print("symbomore4DE yields df+ infos for genes (symbols, biotype...)")
  res.here$ensembl_id <- rownames(res.here)
  # use myobj@row_data 
  rebi <- left_join(as_tibble(res.here) , 
                    myobj@row_data[, c(equi_id, "external_gene_name",
                                       "description", "gene_biotype", 
                                       trustedcolumn)],
                    by=c("ensembl_id"=equi_id) )
  rebi <- rebi %>% relocate(ensembl_id, .before="baseMean")
  rebi[[trustedcolumn]] <- make.unique(rebi[[trustedcolumn]])
  rebi <- dplyr::rename(rebi, gene_name_reps = external_gene_name)
  print(colnames(rebi))
  rebi$gene_name_reps <- unname( sapply(rebi[[trustedcolumn]], function(x) {
    unlist(str_split(x,"\\."))[1]}) )  
  registerduplicates_fromresultsV2(rebi, tabl_dir, outname)
  rebi <- rebi %>% filter(!gene_biotype %in% exclude_types) %>%
    arrange(padj)
  return(rebi)
}

registerduplicates_fromresultsV2 <- function(res_plus_bm, tabl_dir, outname){
  # requires a combined dataframe or tibble DE results + annots (used in symbomore4DE())
  print("SAVing repeated symbols for information: almost always arrives this")
  unique(table(res_plus_bm$gene_name_reps))
  repeatsall <- res_plus_bm %>% group_by(gene_name_reps) %>% 
    summarise(n = n()) %>% filter(n >= 2)
  symbo_bioty <- res_plus_bm[,c("gene_name_reps","gene_biotype")] %>%
    group_by(gene_name_reps) %>%
    dplyr::filter(row_number(gene_name_reps)==1) 
  repeatsall <- left_join(repeatsall, symbo_bioty, by="gene_name_reps")
  write.table(repeatsall, paste0(tabl_dir,"repeatedfeatures_",outname,".tsv"),
              col.names=T, row.names=F, sep='\t' )
  write.table(repeatsall %>% filter(gene_biotype=="protein_coding") ,
              paste0(tabl_dir,"repeatedfeaturesProteinCod_",outname,".tsv"), 
              col.names=T, row.names=F, sep='\t' )
}

joinDEfull_biomart <- function(res.here, bm_sp,equi_id, trustedcolumn,
                               exclude_types, outname, tabl_dir){
  print("if available rnaexp object, use symbomore4DE instead")
  print("joins DE result +biomart annots. Caution:redundant/repeated symbols")
  print("if symbol not found, replaces with ensembl_id")
  res.here$ensembl_id <- rownames(res.here)
  rebi <- left_join(as_tibble(res.here), bm_sp[,c(equi_id, "external_gene_name",
                                           "description", "gene_biotype")],
                    by=c("ensembl_id"=equi_id))
  rebi <- rebi %>% mutate(external_gene_name = case_when(
      external_gene_name == "" ~ ensembl_id,
      is.na(external_gene_name) ~ ensembl_id,
      TRUE ~ external_gene_name ) #end case_when
    )
  rebi <- rebi %>%  relocate(ensembl_id, .before="baseMean")
  rebi[[trustedcolumn]] <- make.unique(rebi$external_gene_name)
  registerduplicates_fromresults(rebi, tabl_dir, outname)
  rebi <- rebi %>% filter(!gene_biotype %in% exclude_types) %>%
    arrange(padj)
  return(rebi)
}


registerduplicates_fromresults <- function(res_plus_bm, tabl_dir, outname){
  # requires a combined dataframe or tibble DE results + biomart annots
  print("SAVing repeated symbols for information: almost always arrives this")
  unique(table(res_plus_bm$external_gene_name))
  repeatsall <- res_plus_bm %>% group_by(external_gene_name) %>% 
    summarise(n = n()) %>% filter(n > 2)
  symbo_bioty <- res_plus_bm[,c("external_gene_name","gene_biotype")] %>%
    group_by(external_gene_name) %>%
    dplyr::filter(row_number(external_gene_name)==1) 
  repeatsall <- left_join(repeatsall, symbo_bioty, by="external_gene_name")
  write.table(repeatsall, paste0(tabl_dir,"repeatedfeatures_",outname,".tsv"),
              col.names=T, row.names=F, sep='\t' )
  write.table(repeatsall %>% filter(gene_biotype=="protein_coding") ,
              paste0(tabl_dir,"repeatedfeaturesProteinCod_",outname,".tsv"), 
              col.names=T, row.names=F, sep='\t' )
}

symbolsasrownames <- function(DEres, trustedcolumn){
  # exclude NA padj or LFC, set new rownames
  # TODO : add error control: trustedcolumn must be IN and unique
  DEres <- DEres %>% filter(!(is.na(padj) | is.na(log2FoldChange)))
  oblidf <-  as.data.frame(DEres[, c(trustedcolumn,"baseMean",  
           "log2FoldChange",   "lfcSE", "stat", "pvalue", "padj" ) ] )
  rownames(oblidf) <- DEres[[trustedcolumn]]
  return(oblidf[,-1])
}

reformat4plot_withlabels <- function(DEres, rtop, sigcutoff,
                                     absLFCcutoff,
                                     trustedcolumn, 
                                 aLFClab, pjlab, col_here){
    genesvizA <- rtop %>% dplyr::filter(abs(log2FoldChange) >= aLFClab  &
                                 padj <= pjlab) %>% pull(trustedcolumn)
  print(genesvizA)
  # slice_min, slice_max : 10 extreme genes (5,5) by default, from rtop
  genesvizB <- rtop %>% slice_min(padj, n=5) %>% pull(trustedcolumn)
  genesvizC <- rtop %>% slice_max(abs(log2FoldChange), n=5) %>% pull(trustedcolumn)
  genesviz <- c(genesvizA, genesvizB, genesvizC)
  DEres$labelrich <- ""
  for (i in genesviz){ DEres[i,]$labelrich <- i  }
  DEres$labelfew <- ""
  for (i in c(genesvizB,genesvizC)){ DEres[i,]$labelfew <- i }
  DEres <- DEres %>% mutate(colorcode = case_when(
    padj <= sigcutoff & log2FoldChange >= absLFCcutoff ~ col_here[3],
    padj <= sigcutoff & log2FoldChange <= -absLFCcutoff ~ col_here[1],
    padj <= sigcutoff & abs(log2FoldChange) < absLFCcutoff ~ col_here[2],
    TRUE ~ "gray30" ) )
  return(DEres)
}

testmycolors <- function(col_here){
  nc <- length(col_here)
  testecolor <- data.frame("a"=seq(nc), 
                           "b"=sapply(seq(nc), function(x) 1/x**2))
  testecolor$col <- col_here
  ggplot(testecolor) + geom_point(aes(a,b,color=col,size=2), alpha=0.6) +
    scale_color_identity()
}

dovolcano_picklabcolumn <- function(rall4plot, labcolumn, strcontrast,
                                    sigcutoff, absLFCcutoff,
                                    pjlab, aLFClab,
                                    outname,  mycaption){
  volca <- ggplot(rall4plot, aes(x=log2FoldChange, y = -log10(padj), 
                                 colour=colorcode)) +
    geom_point( size=.6, alpha=.8) + scale_colour_identity() +
    theme_bw() +
    geom_vline(xintercept = c(absLFCcutoff,-absLFCcutoff),
               color="gray10",
               linetype="dashed", size=.2, alpha=.6) +
    geom_hline(yintercept = -log10(sigcutoff), color= "gray10",
               linetype="dashed", size=.2, alpha=.6) +
    ggrepel::geom_text_repel( aes(label=labcolumn),
                              size=3, max.overlaps=70,  force_pull=5) +
    labs(title=paste(outname, ", ", strcontrast), 
         caption=paste("dashed lines : padj <=", sigcutoff,
                       ", absLFC cutoff :",absLFCcutoff,
                       "\n(",mycaption, ")")
    )
  return(volca)
}

dovolcano_nolabel <- function(rall4plot, strcontrast,
                                    sigcutoff, absLFCcutoff,
                                    pjlab, aLFClab, outname,  mycaption){
  print(outname)
  volca <- ggplot(rall4plot, aes(x=log2FoldChange, y = -log10(padj), 
                                 colour=colorcode)) +
    geom_point( size=.6, alpha=.8) + scale_colour_identity() +
    theme_bw() +
    geom_vline(xintercept = c(absLFCcutoff,-absLFCcutoff),
               color="gray10",
               linetype="dashed", size=.2, alpha=.6) +
    geom_hline(yintercept = -log10(sigcutoff), color= "gray10",
               linetype="dashed", size=.2, alpha=.6) +
    labs(title=paste(outname, ", ", strcontrast), 
         caption=paste("dashed lines : padj <=", sigcutoff,
                       ", absLFC cutoff :",absLFCcutoff,
                       "\n(",mycaption, ")")
    )
  return(volca)
}

givegmt <- function(iwant, myspeci){
  # Example : agmt <- givegmt("CP:REACTOME", "Mus musculus")
  if (iwant == "GO:BP"){
    gmt <- msigdbr(species=myspeci,
                   category='C5',
                   subcategory=iwant)
  }else if (iwant == "CP:KEGG"){
    gmt <- msigdbr(species=myspeci,
                   category='C2',
                   subcategory=iwant)
  }else if (iwant == "CP:REACTOME"){
    gmt <- msigdbr(species=myspeci,
                   category='C2',
                   subcategory=iwant)
  }else{print(iwant);gmt <- NULL}
  return(gmt)
}



paths4plot <- function(fgseaRes, howmanyeach){
  pathup <- fgseaRes[ES>0][head(order(padj),n=howmanyeach),pathway,NES]
  pathdown <- fgseaRes[ES<0][head(order(padj),n=howmanyeach),pathway,NES]
  pathup <- pathup[][order(-NES),pathway]
  pathdown <- pathdown[][order(NES),pathway]
  return(list("down"=pathdown,"up"=pathup))
}

rmduplirows_byensemblidversion <- function(bm_sp){
  # note: invariably highly recommended to use _id_version here
  print("Tipically biomart bm mutiplicates ensembl_gene_id_version  ,")
  print(" yielding df having unique 'ensembl_gene_id_version")
  bm_sp <- bm_sp %>% group_by(ensembl_gene_id_version) %>% 
    dplyr::filter(row_number(ensembl_gene_id_version)==1) 
  return(bm_sp)
}

getorthol_frommart <- function(mart_sp){
  orthoinfo <- getBM(attributes=c( "ensembl_gene_id",
                                    "ensembl_gene_id_version",
                                   "hsapiens_homolog_ensembl_gene",
                                   "hsapiens_homolog_associated_gene_name",
                                   "hsapiens_homolog_orthology_type",
                                   "hsapiens_homolog_orthology_confidence" ),
                     mart = mart_sp) 
  return(orthoinfo)
}

## important: BioMart connect and save into biomart_db locally
savebiomartinfos <- function(SPECIESensmbldsetname, outdir, shortname){
  print("saving biomart frequently used datasets")
  # generates  .rds files, does not overwrites. 
  true_outdir <- paste0(outdir, "biomart_db/")
  if (!dir.exists(true_outdir)){
    dir.create(true_outdir, recursive = FALSE)
  } # end if
  if (file.exists(paste0(true_outdir,"martfull_",shortname, ".rds"))){
    mart_sp <- readRDS(paste0(true_outdir,"martfull_",shortname, ".rds"))
  }else{
    mart_sp <- useEnsembl(biomart = "genes", dataset = SPECIESensmbldsetname)
    saveRDS(mart_sp, paste0(true_outdir,"martfull_",shortname, ".rds"))
  } # end if
  if (file.exists(paste0(true_outdir,"bm_", shortname,".rds"))){
    print("nothing to do")
  } else{  # invariably both _gene_id and _gene_id_version included
    bm_sp <- getBM(attributes=c("ensembl_gene_id",
                                "entrezgene_id",
                                "ensembl_gene_id_version",
                                "external_gene_name",
                                "description", 
                                "gene_biotype"),
                   mart = mart_sp)
    entrezmulti <- bm_sp[,c("entrezgene_id","ensembl_gene_id_version")]
    print("improvement : rmduplirows_byensemblidversion included")
    bm_sp <- rmduplirows_byensemblidversion(bm_sp)
    print("many  ensembl_gene_id_version: RESOLVED (rmduplirow..)")
    saveRDS(bm_sp, paste0(true_outdir,"bm_", shortname,".rds"))
    saveRDS(entrezmulti, paste0(true_outdir,"entrezfull_",shortname,".rds"))
  } # end if
  if (str_to_lower(shortname) != "human"){
    print("getting orthologues list for this no-human species")
    orthofiname <- paste0(true_outdir,"orthoinfo_", shortname,".rds")
    if (file.exists(orthofiname)){
      print("nothing to do, file exists")
    } else{
      orthoinfo <- getorthol_frommart(mart_sp)
      saveRDS(orthoinfo, orthofiname)
    } # end if 
  } #end if str_to_lower
  
  print(paste("END biomart import, and saved into: ", true_outdir))
}

save_ncbi_hits <- function( nosymbolv, shortname){
  # example : save_ncbi_hits( vector_ensembl_id_nosymbol, "rat")
  fiann <- paste0("biomart_db/ncbi_suppl_", shortname, ".tsv")
  if (file.exists(fiann)){
    print(paste("if need a new ncbi list, DELETE : ", fiann))
    Lfil <- opentsv(fiann)
  } else {
    if (shortname == "human"){
      library("org.Hs.eg.db")
      L <- select(org.Hs.eg.db, keys=nosymbolv, keytype="ENSEMBL",
                  columns=c("ENTREZID","SYMBOL", "GENENAME", "GENETYPE"))
    }else if (shortname == "rat"){
      library("org.Rn.eg.db")
      L <- select(org.Rn.eg.db, keys=nosymbolv, keytype="ENSEMBL",
                  columns=c("ENTREZID","SYMBOL", "GENENAME", "GENETYPE"))
    }else { print ("for now, only human and rat ") }
  Lfil <- L %>% filter(!is.na(ENTREZID)) %>% filter(!duplicated(ENSEMBL))
  print("1:many mapping was deduplicated (!duplicated(ENSEMBL))")
  savetsv(Lfil, fiann)
  }
  return(Lfil)
}


selgenes_togg <- function(amat, abm, coluguide, geneslist){
  # input matrix has ensembl_id (not suitable if symbols instead)
  idf <- abm %>% filter(external_gene_name %in% geneslist &
                          gene_biotype == "protein_coding") 
  idf <- idf[,c(coluguide, "external_gene_name")]
  existinginmat <- intersect(rownames(amat), idf[[coluguide]])
  xx <- amat[existinginmat, ]
  #xx[[coluguide]] <- rownames(xx)
  mtdf <- melt(xx, value.name = "expression_values")

  colnames(mtdf) <- c(coluguide, "sample","expression_values")
  mtdf <- left_join(mtdf, idf)
  mygg <- ggplot(mtdf, aes(sample,expression_values, 
                           color=external_gene_name) ) +
    geom_point() +   
    scale_color_brewer(type="qual",palette=7) +
    theme_bw() + coord_flip() +
    labs(title="Selected genes: expression across samples")  
  return(mygg)
}

doboxplot3PC_bicond <- function(pca.res, metadata, TWOCOLORS,mytitle){
  metadata$condition <- as.factor(metadata$condition)
  pic <- data.frame(pca.res$ind$coord)
  pic$sample <- rownames(pca.res$ind$coord)
  pic <- inner_join(pic, metadata, by="sample")
  dimsme <- melt(pic[,c("Dim.1","Dim.2","Dim.3","condition")])
  agg <- ggplot(dimsme, aes(x=condition, y=value, col=condition)) +
    geom_point(position=position_jitterdodge()) + 
    scale_colour_manual(values=TWOCOLORS) +
    geom_boxplot(alpha=0.7) + theme_light() + facet_grid(~variable) +
    labs(title=mytitle)
  return(agg)
}


plotPCAvarswithsign <- function(res.pca, mypval, metadata,
                                CHOSENCOLORS, tooextreme=c(""), invcol=FALSE){
  coeffi <- dimdesc(res.pca, axes=1)
  coefdf <- as.data.frame(coeffi$Dim.1$quanti)
  coefdf$symbol_unique <- rownames(coefdf)
  print(tooextreme)
  exfiltp <- coefdf %>% filter(!symbol_unique %in% tooextreme) %>%
    filter(p.value <= mypval)
  if(invcol){
    colarrow <- rev(CHOSENCOLORS)
  }else{colarrow <- CHOSENCOLORS}
  colorcorr <- ifelse(res.pca$var$coord[,"Dim.1"] > 0, colarrow[1],colarrow[2])
  colorcorr <- factor(colorcorr, levels=colarrow )
  theggpl <- fviz_pca_biplot(res.pca, # check
                             geom.ind="point",
                             fill.ind=metadata$condition, #check
                             col.ind="white",
                             pointshape=21,
                             pointsize=3,
                             labelsize=2,
                             palette=CHOSENCOLORS,  #check
                             addEllipses=FALSE,
                             invisible="quali",
                             repel=TRUE,
                             col.var=colorcorr,
                             alpha.var="contrib",
                             select.var=list(names=exfiltp$symbol_unique)) +
    labs(fill="Condition",color="Sign.Dim1",alpha="Contrib", caption=mypval) +
    scale_color_manual(name = "Sign.Dim1", labels = c("+","-"),
                       values=rev(colarrow))
  print("function about to end")
  return(theggpl)
}

printexampleggplottext <- function(){
  mye = "ggplot() + annotate('text', x=1,y=1,size = 4,label = 'MYBLAH') + theme_void()"
  print(mye)
}

getcoefvar <- function(amat){
  cat("getcoefvar(amat) outputs a dataframe with cols: $symbol_unique $coef.var\n")
  ameans <-apply(amat, 1, mean)
  avars <- apply(amat, 1, var)
  coefv <- sqrt(avars) / (ameans+1)
  coefv[is.na(coefv)] <- 0
  return(tibble("symbol_unique"=rownames(amat),"coef.var"=coefv))
} # ok, ids replaced by symbol_unique does not change anything :)

coefvar_bycond <- function(amat, ametadata){
  cat("coefvar_bycond(amat,ametadata) outputs dataframe: n conditions\n")
  condis <- unique(ametadata$condition)
  grs <- list(); cvdfs <- data.frame("symbol_unique"=rownames(amat))
  for (i in condis){
    grs[[i]] <- ametadata %>% filter(condition==i) %>% pull(sample)
  }
  for (i in condis){
    submatr <- amat[, grs[[i]]]
    icoefdf <- getcoefvar(submatr)
    colnames(icoefdf) <- c("symbol_unique",i)
    cvdfs <- merge(cvdfs,icoefdf, by="symbol_unique", sort=FALSE)
  }
  return(cvdfs)
}

getdivergent4complexheat <- function(logFCvector){
  cat(" getdivergent4complexheat :  outputs a function")
  cat("   (needs ComplexHeatmap and circlize)\n")
  n = length(logFCvector)
  maxabs = max(abs(logFCvector))
  varifun = circlize::colorRamp2(seq(-maxabs, maxabs, length = n), 
                                 hcl.colors(n,"Blue-Red 2"))
  return(varifun)
}

geometricMean_bycond <- function(amat, ametadata, eps){
  cat('geometricMean_bycond(amat, ametadata) outs df: symbol_unique, n conditions')
  cat('if only one condition, use just getgeommean \n')
  minvalnonzero <- min(amat[amat > 0])
  condis <- unique(ametadata$condition)
  grs <- list(); gmdfs <- data.frame("symbol_unique"=rownames(amat))
  for (i in condis){
    grs[[i]] <- ametadata %>% filter(condition==i) %>% pull(sample)
  }
  for (i in condis){
    submatr <- amat[,grs[[i]]]
    igmdf <- getgeommean(submatr, minvalnonzero, eps)
    colnames(igmdf) <- c("symbol_unique", i)
    gmdfs <- merge(gmdfs, igmdf , by="symbol_unique", sort=FALSE)
  } 
  return(gmdfs)
}
