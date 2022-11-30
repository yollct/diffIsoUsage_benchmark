# copyright: Xi Wang (xi.wang@newcastle.edu.au)

######################
## Class SeqGeneSet ##
######################

SeqGeneSet <- setClass("SeqGeneSet", 
                       representation = representation(
                         name = "character", 
                         sourceFile = "character", 
                         geneList = "character", 
                         GS = "list",  
                         GSNames = "character", 
                         GSDescs = "character", 
                         GSSize = "numeric", 
                         GSSizeMin = "numeric",
                         GSSizeMax = "numeric", 
                         GS.Excluded = "list",
                         GSNames.Excluded = "character",
                         GSDescs.Excluded = "character", 
                         GSEA.ES = "numeric",
                         GSEA.ES.pos = "numeric",
                         GSEA.ES.perm = "matrix",
                         GSEA.score.cumsum = "matrix", 
                         GSEA.normFlag = "logical", 
                         GSEA.pval = "numeric",
                         GSEA.FWER = "numeric", 
                         GSEA.FDR = "numeric", 
                         sc.ES = "matrix", 
                         sc.ES.perm = "matrix", 
                         sc.normFlag = "logical", 
                         scGSEA = "logical", 
                         sc.pval = "matrix", 
                         sc.FWER = "matrix", 
                         sc.FDR = "matrix", 
                         version = "Versions"), 
                       prototype = prototype(
                         name = new("ScalarCharacter", NA),
                         version = new("Versions", "0.0.2"))
)

newGeneSets <- function(GS, GSNames, GSDescs, geneList, scGSEA = FALSE, 
                        name = NA_character_, sourceFile = NA_character_, 
                        GSSizeMin = 5, GSSizeMax = 1000) {
  stopifnot( length(GS) == length(GSNames) & length(GS) == length(GSDescs)) 
  stopifnot( max(unlist(GS)) <= length(geneList) )
  GSSize = unlist(lapply(GS, length))
  GSUse = GSSize >= GSSizeMin & GSSize <= GSSizeMax
  
  GS <- new("SeqGeneSet", name = name, sourceFile = sourceFile, geneList = geneList, 
            GS = GS[GSUse], GSNames = GSNames[GSUse], GSDescs = GSDescs[GSUse], 
            GSSizeMin = GSSizeMin, GSSizeMax = GSSizeMax, GSSize = GSSize[GSUse], 
            GS.Excluded = GS[!GSUse], GSNames.Excluded = GSNames[!GSUse], GSDescs.Excluded = GSDescs[!GSUse], 
            GSEA.ES = numeric(0), GSEA.ES.pos = numeric(0), 
            GSEA.ES.perm = matrix(0, 0 ,0),  GSEA.score.cumsum = matrix(0, 0 ,0), 
            GSEA.normFlag = FALSE,  GSEA.pval = numeric(0), 
            GSEA.FWER = numeric(0), GSEA.FDR = numeric(0),
            sc.ES = matrix(0, 0 ,0), sc.ES.perm = matrix(0, 0 ,0),
            sc.normFlag = FALSE, scGSEA = scGSEA, 
            sc.pval = matrix(0, 0 ,0), sc.FWER = matrix(0, 0 ,0), sc.FDR = matrix(0, 0 ,0) )
  GS
}

setMethod("show", 
          signature(object="SeqGeneSet"),
          function(object) {
            if(!object@scGSEA) {
              expr1 <- ifelse(length(object@GSEA.ES) > 0, 
                            paste(selectSome(object@GSEA.ES, maxToShow=4), collapse=", "), 
                            "not computed")
              expr2 <- ifelse(length(object@GSEA.ES.pos) > 0, 
                            paste(selectSome(object@GSEA.ES.pos, maxToShow=4), collapse=", "), 
                            "not computed")
              expr3 <- ifelse(ncol(object@GSEA.ES.perm) > 0, 
                            paste(ncol(object@GSEA.ES.perm), "-time permutation", sep=""),  
                            "not performed")
              expr4 <- ifelse(object@GSEA.normFlag, "Yes", "No")
              expr5 <- ifelse(length(object@GSEA.pval) > 0, 
                            paste(selectSome(object@GSEA.pval, maxToShow=4), collapse=", "), 
                            "not computed")
              expr6 <- ifelse(length(object@GSEA.FWER) > 0, 
                            paste(selectSome(object@GSEA.FWER, maxToShow=4), collapse=", "), 
                            "not computed")
              expr7 <- ifelse(length(object@GSEA.FDR) > 0, 
                            paste(selectSome(object@GSEA.FDR, maxToShow=4), collapse=", "), 
                            "not computed")
              cat("SeqGeneSet object: ", object@name, "\n",
                "GeneSetSourceFile: ", object@sourceFile, "\n", 
                "GeneSets: ", paste(selectSome(object@GSNames, maxToShow=4), 
                                    collapse="\n          "), "\n", 
                "  with the number of genes in respective sets: ", 
                paste(selectSome(object@GSSize, maxToShow=4), collapse=", "), "\n",
                "  brief descriptions: \n          ",
                paste(selectSome(object@GSDescs, maxToShow=4), collapse="\n          "), "\n",
                "  # gene sets passed filter: ", size(object), 
                "  (#genes >= ", object@GSSizeMin, " AND <= ", object@GSSizeMax, ")\n", 
                "  # gene sets excluded: ", sizeOfExcludedGS(object), 
                "  (#genes < ", object@GSSizeMin, " OR > ", object@GSSizeMax, ")\n",
                "ES scores: ", expr1, "\n",
                "ES postions: ", expr2, "\n",
                "Permutated ES scores: ", expr3, "\n", 
                "ES scores normalized: ", expr4, "\n",
                "ES p-value: ", expr5, "\n",
                "ES FWER: ", expr6, "\n",
                "ES FDR: ", expr7, "\n",
                sep="")
            } else {
              expr1 <- ifelse(length(object@sc.ES) > 0, 
                              paste(nrow(object@sc.ES), "genesets by", ncol(object@sc.ES), "cells computed"), 
                              "not computed")
              expr2 <- ifelse(length(object@sc.ES.perm) > 0, 
                              paste(nrow(object@sc.ES), "genesets by", ncol(object@sc.ES), "permutation times computed"), 
                              "not computed")
              expr3 <- ifelse(object@sc.normFlag, "Yes", "No")
              expr4 <- ifelse(length(object@sc.pval) > 0, 
                              paste(nrow(object@sc.pval), "genesets by ", ncol(object@sc.pval), " cells computed"), 
                              "not computed")
              expr5 <- ifelse(length(object@sc.FWER) > 0, 
                              paste(nrow(object@sc.FWER), "genesets by ", ncol(object@sc.FWER), " cells computed"), 
                              "not computed")
              expr6 <- ifelse(length(object@sc.FDR) > 0, 
                              paste(nrow(object@sc.FDR), "genesets by ", ncol(object@sc.FDR), " cells computed"), 
                              "not computed")
              
              cat("SeqGeneSet object: ", object@name, "\n",
                "GeneSetSourceFile: ", object@sourceFile, "\n", 
                "GeneSets: ", paste(selectSome(object@GSNames, maxToShow=4), 
                                    collapse="\n          "), "\n", 
                "  with the number of genes in respective sets: ", 
                paste(selectSome(object@GSSize, maxToShow=4), collapse=", "), "\n",
                "  brief descriptions: \n          ",
                paste(selectSome(object@GSDescs, maxToShow=4), collapse="\n          "), "\n",
                "  # gene sets passed filter: ", size(object), 
                "  (#genes >= ", object@GSSizeMin, " AND <= ", object@GSSizeMax, ")\n", 
                "  # gene sets excluded: ", sizeOfExcludedGS(object), 
                "  (#genes < ", object@GSSizeMin, " OR > ", object@GSSizeMax, ")\n",
                "sc-ES scores: ", expr1, "\n",
                "Permutated sc-ES scores: ", expr2, "\n", 
                "sc-ES scores normalized: ", expr3, "\n",
                "sc-ES p-value: ", expr4, "\n",
                "sc-ES FWER: ", expr5, "\n",
                "sc-ES FDR: ", expr6, "\n",
                sep="")
            }
          })

setMethod("[", 
          signature = signature(
            x="SeqGeneSet", i="numeric"),
          function(x, i, j, ..., drop=TRUE) {
            stopifnot(all(i <= size(x)))
            if (length(i) == 1 && drop) 
              return(as.vector(x@GS[[i]]))
            GS <- x
            GS@GS <- x@GS[i]
            GS@GSNames <- x@GSNames[i]
            GS@GSDescs <- x@GSDescs[i]
            GS@GSSize <- x@GSSize[i]
            GS@GS.Excluded <- list()
            GSNames.Excluded <- character()
            GSDescs.Excluded <- character()
            if(length(x@GSEA.ES) > 0 ) 
              GS@GSEA.ES <- x@GSEA.ES[i]
            if(length(x@GSEA.ES.pos) > 0)
              GS@GSEA.ES.pos <- x@GSEA.ES.pos[i]
            if(nrow(x@GSEA.ES.perm) > 0)
              GS@GSEA.ES.perm <- x@GSEA.ES.perm[i,,drop=FALSE]
            if(nrow(x@GSEA.score.cumsum) > 0)
              GS@GSEA.score.cumsum <- x@GSEA.score.cumsum[i,,drop=FALSE]
            if(length(x@GSEA.pval) > 0)
              GS@GSEA.pval <- x@GSEA.pval[i]
            if(length(x@GSEA.FDR) > 0)
              GS@GSEA.FDR <- x@GSEA.FDR[i]
            if(length(x@GSEA.FWER) > 0)
              GS@GSEA.FWER <- x@GSEA.FWER[i]
            if(nrow(x@sc.ES) > 0)
              GS@sc.ES <- x@sc.ES[i,,drop=FALSE]
            if(nrow(x@sc.ES.perm) > 0)
              GS@sc.ES.perm <- x@sc.ES.perm[i,,drop=FALSE]
            if(nrow(x@sc.pval) > 0)
              GS@sc.pval <- x@sc.pval[i,,drop=FALSE]
            if(nrow(x@sc.FWER) > 0)
              GS@sc.FWER <- x@sc.FWER[i,,drop=FALSE]
            if(nrow(x@sc.FDR) > 0)
              GS@sc.FDR <- x@sc.FDR[i,,drop=FALSE]
            GS
          } )

geneSetNames <- function(GS) {
  stopifnot( is(GS, "SeqGeneSet"))
  gsn <- GS@GSNames
  names(gsn) <- GS@GSNames
  gsn
}
  
geneSetDescs <- function(GS) {
  stopifnot( is(GS, "SeqGeneSet"))
  gsd <- GS@GSDescs
  names(gsd) <- GS@GSNames
  gsd
}

geneSetSize <- function(GS) {
  stopifnot( is(GS, "SeqGeneSet"))
  gss <- GS@GSSize
  names(gss) <- GS@GSNames
  gss
}

size <- function(GS) {
  stopifnot( is(GS, "SeqGeneSet"))
  length(GS@GS)
}

geneList <- function(GS) {
  stopifnot( is(GS, "SeqGeneSet"))
  GS@geneList
}


sizeOfExcludedGS <- function(GS) {
  stopifnot( is (GS, "SeqGeneSet") )
  length(GS@GS.Excluded)
}

########################
## Class ReadCountSet ##
########################

#require(Biobase)
ReadCountSet <- setClass( "ReadCountSet", 
          contains = "eSet",
          representation = representation(
            featureData_gene = "data.frame",
            permute_NBstat_exon = "matrix",
            permute_NBstat_gene = "matrix"
            ),
          prototype = prototype( new( "VersionedBiobase",
                                      versions = c( classVersion("eSet"), 
                                                    ReadCountSet = "1.0.0" ) ) )
)

newReadCountSet <- function( readCounts, exonIDs, geneIDs ) {
  readCounts <- as.matrix( readCounts )
  if( any( round( readCounts ) != readCounts ) ){
    stop( "The count data is not integer." )}
  mode( readCounts ) <- "integer"
  nRow <- nrow(readCounts)
  
  phenoData <- annotatedDataFrameFrom( readCounts, byrow=FALSE )
  featureData <- annotatedDataFrameFrom( readCounts, byrow=TRUE )
  
  samples <- colnames(readCounts)
  phenoData$label <- rep(1, length(samples))
  phenoData$label[regexpr("C", samples)==1] <- 0
  phenoData$label <- as.factor(phenoData$label)
  varMetadata( phenoData )[ "label", "labelDescription" ] <- "label of samples (1's for cases, 0's for controls" 
  
  exonIDs <- as.character( exonIDs )
  if( length(exonIDs) != nRow )
    stop( "exonIDs must be of the same length as the number of columns in readCounts" )
  geneIDs <- as.factor( geneIDs )
  if( length(geneIDs) != nRow )
    stop( "geneIDs must be of the same length as the number of columns in readCounts" )  
  featureData$exonIDs <- exonIDs
  varMetadata( featureData )[ "exonIDs", "labelDescription" ] <- "exon ID (unique only within a gene)" 
  featureData$geneIDs <- geneIDs
  varMetadata( featureData )[ "geneIDs", "labelDescription" ] <- "ID of gene to which the exon belongs"
  
  featureData$testable <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "testable", "labelDescription" ] <- "slot indicating if an exon should be considered in the test"
  featureData$prob_case <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "prob_case", "labelDescription" ] <- "expection of this exon's prob in a gene in the case group"
  featureData$prob_ctrl <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "prob_ctrl", "labelDescription" ] <- "expection of this exon's prob in a gene in the control group"
  featureData$var_case <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "var_case", "labelDescription" ] <- "variance of this exon's pron in the case group"
  featureData$var_ctrl <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "var_ctrl", "labelDescription" ] <- "variance of this exon's pron in the control group"
  featureData$NBstat <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "NBstat", "labelDescription" ] <- "NB stat for DE"
  featureData$pvalue <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "pvalue", "labelDescription" ] <- "p-value from NB permutation"
  featureData$padjust <- rep( NA_real_, nRow )
  varMetadata( featureData )[ "padjust", "labelDescription" ] <- "adjusted p-value from NB permutation"
  
  featureData_gene_NBstat <- rep( NA_real_, length(unique(geneIDs)))
  featureData_gene_pval <- rep( NA_real_, length(unique(geneIDs))) 
  featureData_gene_padj <- rep( NA_real_, length(unique(geneIDs))) 
    
  RCS <- new( "ReadCountSet",
              assayData = assayDataNew( "environment", counts=readCounts ),
              phenoData = phenoData,
              featureData = featureData,
              featureData_gene = data.frame(NBstat = featureData_gene_NBstat,
                                            pval = featureData_gene_pval,
                                            padj = featureData_gene_padj,
                                            row.names = unique(as.character(geneIDs))),
              permute_NBstat_exon = matrix(0,0,0),
              permute_NBstat_gene = matrix(0,0,0)
              )
  RCS
}

setValidity( "ReadCountSet", function( object ) {  
 if( ! "geneIDs" %in% names(fData(object)) )
    return( "featureData does not contain a 'geneIDs' column.")
  if( ! is( fData(object)$geneID, "factor" ) )
    return( "The 'geneID' column in fData is not a factor." )
  if( ! "exonIDs" %in% names(fData(object)) )
    return( "featureData does not contain an 'exonIDs' column.")
  if( ! is( fData(object)$exonIDs, "character" ) )
    return( "The 'exonIDs' column in fData is not a character vector." )
 if( !is(fData(object)$NBstat, "numeric")){
   return( "The 'NBstat' values are not numeric")}
 if( !is(fData(object)$pvalue, "numeric")){
   return( "The 'pvalue' values are not numeric")}
 if( !is(fData(object)$padjust, "numeric")){
   return( "The 'padjust' values are not numeric")} 
 if( !is.integer( assayData(object)[["counts"]] ) )
   return( "The count data is not in integer mode." )
 if( any( assayData(object)[["counts"]] < 0 ) )
   return( "The count data contains negative values." )
 if( !is(object@featureData_gene , "data.frame")){
   return( "The 'featureData_gene' is not data.frame")}
 if( !is(object@permute_NBstat_exon , "matrix")){
   return( "The 'permute_NBstat_exon' is not matrix")}
 if( !is(object@permute_NBstat_gene , "matrix")){
   return( "The 'permute_NBstat_gene' is not matrix")}
 
 TRUE
}) 

loadExonCountData <- function(case.files, control.files) {
  nCase <- length(case.files)
  nCtrl <- length(control.files)
  sample.names <- c(paste('S',1:nCase,sep=""), paste('C',1:nCtrl,sep=""))
  lf <- lapply(c(case.files, control.files), function(x) read.table(x, header = FALSE, 
                                                                    stringsAsFactors = FALSE))
  if (!all(sapply(lf[-1], function(x) all(x$V1 == lf[1]$V1)))) 
    stop("Count files have differing gene ID column.")
  dcounts <- sapply(lf, `[[`, "V2")
  rownames(dcounts) <- lf[[1]][, 1]
  colnames(dcounts) <- sample.names
  
  # dcounts for EXONs, and gcounts for GENEs
  dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] 
  geneIDs <- sapply(strsplit(rownames(dcounts), ":"), "[[", 1)
  exonIDs <- paste("E", sapply(strsplit(rownames(dcounts), ":"), function(x) {
                               return(x[2])}), sep = "")
  #gcounts <- do.call(rbind, tapply(1:nrow(dcounts), as.factor(geneIDs), function(rows) {colSums(dcounts[rows, , drop = FALSE]) } ))
    print(head(dcounts))
  idx <- order(geneIDs) # to make sure gene IDs are sorted in lexicographic order; otherwise will cause error in applying `tapply`
  RCS <- newReadCountSet( dcounts[idx, ], exonIDs[idx], geneIDs[idx] )
  RCS 
}

label <- function( RCS ) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  lb <- pData(RCS)$label
  names( lb ) <- colnames( counts(RCS) )
  lb
}

setMethod("counts", signature(object="ReadCountSet"), 
  function( object ) {
    RCS <- object
    assayData(RCS)[["counts"]]
})

setReplaceMethod("counts", signature(object="ReadCountSet", value="matrix"),
  function( object, value ) {
    RCS <- object
    assayData(RCS)[[ "counts" ]] <- value
    validObject(RCS)
    RCS
})

geneID <- function( RCS ) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  g <- fData(RCS)$geneIDs
  names(g) <- rownames( counts(RCS) )
  g
}

`geneID<-` <- function( RCS, value ) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  fData(RCS)$geneIDs <- as.factor(value)
  validObject(RCS)
  RCS
}

exonID <- function( RCS ) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  g <- fData(RCS)$exonIDs
  names(g) <- rownames( counts(RCS) )
  g
}

`exonID<-` <- function( RCS, value ) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  fData(RCS)$exonIDs <- as.character(value)
  validObject(RCS)
  RCS
}

subsetByGenes <- function( RCS, genes ) {
  stopifnot( is( RCS, "ReadCountSet" ) )
  stopifnot( all( genes %in% levels(geneID(RCS)) ) )
  indxExon <- as.character(geneID(RCS)) %in% genes
  indxGene <- unique(as.character(geneID(RCS))) %in% genes
  RCS2 <- RCS[ indxExon, ]
  fData(RCS2)$geneIDs <- as.factor(as.character(fData(RCS2)$geneIDs))
  RCS2@featureData_gene <- RCS@featureData_gene[indxGene,]
  if( length(RCS@permute_NBstat_exon) ) {
    RCS2@permute_NBstat_exon <- as.matrix(RCS@permute_NBstat_exon[indxExon, ])
    RCS2@permute_NBstat_gene <- as.matrix(RCS@permute_NBstat_gene[indxGene, ])
  } else {
    RCS2@permute_NBstat_exon <- matrix(0,0,0)
    RCS2@permute_NBstat_gene <- matrix(0,0,0)
  }
  RCS2
}