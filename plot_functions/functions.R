

colMap <- c(iso_ktsp="#999999", drimseq="#E69F00", dexseq="#56B4E9", dturtle="#783f04", seqGSEA="#F0E442", cuffdiff="#0072B2", junctionseq="#D55E00", saturn="#CC79A7", DSGseq="#9ACD32", nbsplice="#00868B", LimmaDS="#5569da", edgeR="#B536DA")

list_false_positive <- function(methoddf, truthfile, method){
    met <- methoddf[,method]
    names(met) <- methoddf$feature_id
    fp <- names(met)[met<0.05][!names(met)[met<0.05] %in% truthfile[truthfile$status==1,]$feature_id]
    return(fp)
}

list_false_negative <- function(methoddf, truthfile, method){
    met <- methoddf[,method]
    names(met) <- methoddf$feature_id
    fn <- truthfile[truthfile$status==1,]$feature_id[truthfile[truthfile$status==1,]$feature_id %in% names(met)[!met<0.05]]
    return(fn)
}

list_tp <- function(methoddf, truthfile, method){
    met <- methoddf[,method]
    names(met) <- methoddf$feature_id
    tp <- unique(names(met)[met<0.05])[unique(names(met)[met<0.05]) %in% unique(truthfile[truthfile$status==1,]$feature_id)]
    return(tp)
}

# list_tp(methoddf_g, truthfile_g, "dexseq_stageR") %>% length
# list_false_positive(methoddf_g, truthfile_g, "dexseq_stageR") 
# list_false_negative(methoddf_g, truthfile_g, "dexseq_stageR")

# tovenn = list(dexseq= names(met)[met<0.05],
#             truth = truthfile[truthfile$status==1,]$feature_id)
# library(ggvenn)
# ggvenn(tovenn)

#calculate fdr
cal_fdr_tpr <- function(methoddf, truthfile, thresholds){
    mets <- colnames(methoddf)[2:ncol(methoddf)]
    out <- do.call(cbind, sapply(mets, function(method){
        met <- methoddf[,method]
        names(met) <- methoddf$feature_id
        fdrs <- sapply(thresholds, function(x){
            fp <- sum(!names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
            pp <- length(names(met)[met<x])
            fp/pp})

        
        tprs <- sapply(thresholds, function(x){
            tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
            fn <- sum(truthfile[truthfile$status==1,]$feature_id %in% names(met)[met>=x])
            p <- length(truthfile[truthfile$status==1,]$feature_id)
            tp/(tp+fn)})

        outout <- data.frame(fdr=c(fdrs, tprs), type=rep(c("fdr", "tpr"), each=length(fdrs)), tool=rep(method, length(fdrs)*2), thresholds=rep(thresholds, 2))
        list(t(outout))
        }
    ))
    

    return(out)
}


#calculate precision recall
cal_pre_re <- function(methoddf, truthfile, split=NULL, tx=TRUE){
    if (is.null(split)){
        mets <- colnames(methoddf)[2:ncol(methoddf)]
        xx<-1
        out <- do.call(cbind, sapply(mets, function(method){
            methoddf <- methoddf[methoddf$feature_id %in% truthfile$feature_id,]
            namet <- methoddf[,method]
            
            names(namet) <- "value"
            namet$feature_id <- methoddf$feature_id
            x<-thresholds[method,]$thres
            

            namet$value <- lapply(namet$value, function(y){ifelse(x==0.05, ifelse(is.na(y), 1, y), ifelse(is.na(y), 0, y))}) %>% unlist()
            
            #precisions <- sapply(thresholds, function(x){
            boomet <- lapply(namet$value, function(y){ifelse(x==0.05, y<x, y>x)}) %>% unlist

            tp <- sum(unique(namet$feature_id[boomet]) %in% truthfile[truthfile$status==1,]$feature_id) 
            fp <- sum(!unique(namet$feature_id[boomet]) %in% truthfile[truthfile$status==1,]$feature_id) 
            fn <- sum(!unique(namet$feature_id[!boomet]) %in% truthfile[truthfile$status==0,]$feature_id) 
        
            pp <- length(unique(namet$feature_id[boomet]))
            precisions <- tp/pp

             #   })
            #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
            #recall <- sapply(thresholds, function(x){
            
            p <- length(truthfile[truthfile$status==1,]$feature_id)
            recall <- tp/p

            f1 <- 2*tp/(2*tp+fp+fn)
            #    return(tp/p)})
            xx <- xx+1
            outout <- data.frame(fdr=c(precisions, recall,f1, tp, p, pp), type=rep(c("precision", "recall","f1", "truepos", "realpos", "detectedpos"), each=length(precisions)), tool=rep(method, length(precisions)*3), thresholds=rep(x, 3))
            list(t(outout))
            }
        ))
    } else {
        out <- data.frame(row.names=c("fdr","type","tool","thresholds","splits"))
        if (split=='events'){
            splitcat <- c("DTE", "DTU", "IS")
        } else {
            splitcat <- unique(truthfile[split][,1])
        }
        print(splitcat)
        for (i in splitcat){
            splitted_truth <- truthfile %>% dplyr::filter(truthfile[split]==i)
            
            not_this_truth <- truthfile[truthfile[split]!=i,]$feature_id 

            
            ## filter gene that is not this category
            splitted_method <- methoddf %>% dplyr::filter(!feature_id %in% not_this_truth)
            
            
            thisout <- do.call(cbind, sapply(colnames(splitted_method)[2:ncol(splitted_method)], function(method){
                namet <- splitted_method[,method]
                names(namet) <- "value"
                namet$feature_id <- splitted_method$feature_id
                
                x <- thresholds[method,]$thres
                namet$value <- lapply(namet$value, function(y){ifelse(x==0.05, ifelse(is.na(y), 1, y), ifelse(is.na(y), 0, y))}) %>% unlist()
                print(method)
                boomet <- lapply(namet$value, function(y){ifelse(x==0.05, y<x, y>x)}) %>% unlist
                #precisions <- sapply(thresholds, function(x){
                
                tp <- sum(unique(namet$feature_id[boomet]) %in% splitted_truth[splitted_truth$status==1,]$feature_id) 
                fp <- sum(!unique(namet$feature_id[boomet]) %in% splitted_truth[splitted_truth$status==1,]$feature_id) 
                fn <- sum(!unique(namet$feature_id[!boomet]) %in% splitted_truth[splitted_truth$status==0,]$feature_id) 
                
                pp <- length(unique(namet$feature_id[boomet]))
                print(pp)
                precisions <- tp/pp
                f1 <- 2*tp/(2*tp+fp+fn)
                #})
                #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
                #recall <- sapply(thresholds, function(x){
                if (tx){
                    tp <- sum(namet$feature_id[boomet] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                    p <- length(splitted_truth[splitted_truth$status==1,]$feature_id)
                } else {
                    tp <- sum(namet$feature_id[boomet] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                    p <- length(splitted_truth[splitted_truth$status==1,]$feature_id)
                }
            
                
                recall <- tp/p
    
            #})
                outout <- data.frame(fdr=c(precisions, recall, f1, tp, p, pp), type=rep(c("precision", "recall", "f1", "truepos", "realpos", "detectedpos"), each=length(precisions)), tool=rep(method, length(precisions)*3), thresholds=rep(x,3), splits=rep(i, length(precisions)*3))
                list(t(outout))
            }))
            out <- cbind(out, thisout)
        }
    }

    return(out)
}

cal_f1 <- function(methoddf, truthfile, split=NULL, tx=TRUE){
    if (is.null(split)){
        mets <- colnames(methoddf)[2:ncol(methoddf)]
        xx<-1
        out <- do.call(cbind, sapply(mets, function(method){
            namet <- methoddf[,method]
            names(namet) <- methoddf$feature_id
            x<-thresholds[method,]$thres
          
            met <- lapply(namet, function(y){ifelse(x==0.05, ifelse(is.na(y), 1, y), ifelse(is.na(y), 0, y))}) %>% unlist()
            #precisions <- sapply(thresholds, function(x){
            boomet <- lapply(met, function(y){ifelse(x==0.05, y<x, y>x)}) %>% unlist
            
            tp <- sum(names(met)[boomet] %in% truthfile[truthfile$status==1,]$feature_id) 
            fp <- sum(!names(met)[boomet] %in% truthfile[truthfile$status==1,]$feature_id) 
            fn <- sum(!names(met)[!boomet] %in% truthfile[truthfile$status==0,]$feature_id) 
            
            pp <- length(names(met)[boomet])
            f1 <- 2*tp/(2*tp+fp+fn)
             #   })
            #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
            #recall <- sapply(thresholds, function(x){
            
            #    return(tp/p)})
            xx <- xx+1
            outout <- data.frame(f=f1, type="f1", tool=method)
            list(t(outout))
            }
        ))
    } else {
        out <- data.frame(row.names=c("fdr","type","tool","thresholds","splits"))
        for (i in unique(truthfile[split][,1])){
            splitted_truth <- truthfile %>% dplyr::filter(truthfile[split]==i)
            not_this_truth <- truthfile[truthfile[split]!=i,]$feature_id    
            splitted_method <- methoddf %>% dplyr::filter(!feature_id %in% not_this_truth)
            thisout <- do.call(cbind, sapply(colnames(splitted_method)[2:ncol(splitted_method)], function(method){
                namet <- splitted_method[,method]
                names(namet) <- splitted_method$feature_id
                x <- thresholds[method,]$thres
                met <- lapply(namet, function(y){ifelse(x==0.05, ifelse(is.na(y), 1, y), ifelse(is.na(y), 0, y))}) %>% unlist()

                boomet <- lapply(met, function(y){ifelse(x==0.05, y<x, y>x)}) %>% unlist
                #precisions <- sapply(thresholds, function(x){
                tp <- sum(names(met)[boomet] %in% truthfile[truthfile$status==1,]$feature_id) 
                fp <- sum(!names(met)[boomet] %in% truthfile[truthfile$status==1,]$feature_id) 
                fn <- sum(!names(met)[!boomet] %in% truthfile[truthfile$status==0,]$feature_id) 
                print(tp)
                pp <- length(names(met)[boomet])
                f1 <- 2*tp/(2*tp+fp+fn)
                #})
                #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
                #recall <- sapply(thresholds, function(x){
            
            
            #})
                outout <- data.frame(fdr=f1, type="f1", tool=method, splits=split)
                list(t(outout))
            }))
            out <- cbind(out, thisout)
        }
    }

    return(out)
}


plot_upset <- function(methoddf, thresholds){
    x<-1
    metname <- colnames(methoddf)[2:ncol(methoddf)]
    thisout <- do.call(cbind, sapply(metname, function(method){
        met <- methoddf[method]
        met$new <- ifelse(met<thresholds, 1, 0)
        x<-x+1
        data.frame(met$new)
    }))
    thisout
}

pivot_output <- function(outputpr, split=NULL){
    outputpr <- as.data.frame(t(outputpr))
    wideoutputpr <- pivot_wider(outputpr, names_from="type", values_from="fdr")
    wideoutputpr$precision <- as.numeric(wideoutputpr$precision)
    wideoutputpr$recall <- as.numeric(wideoutputpr$recall)
    if (is.null(split)) {
        wideoutputpr$tool <- factor(wideoutputpr$tool)
    } else if (split=="gene_group") {
        wideoutputpr <- wideoutputpr %>% dplyr::filter(splits != "1")
        wideoutputpr$tool <- factor(wideoutputpr$tool)
        wideoutputpr$splits <- factor(wideoutputpr$splits, level=c("2-4","5-9", ">9")) 
        wideoutputpr$recall[is.nan(wideoutputpr$recall)] = 0
    } else if (split=="fc_group") {
        wideoutputpr <- wideoutputpr %>% dplyr::filter(splits != "1")
        wideoutputpr$tool <- factor(wideoutputpr$tool)
        wideoutputpr$splits <- factor(wideoutputpr$splits, level=c("2","3","4","5"))
        wideoutputpr$recall[is.nan(wideoutputpr$recall)] = 0
    } else if (split=="events") {
        wideoutputpr <- wideoutputpr %>% dplyr::filter(splits != "0")
        wideoutputpr$tool <- factor(wideoutputpr$tool)
        wideoutputpr$splits <- factor(wideoutputpr$splits, level=c("DTE","DTU","IS"))
        wideoutputpr$recall[is.nan(wideoutputpr$recall)] = 0
    } 

    return(wideoutputpr)
}

compare_tools <- function(x,y,z=NULL,split=FALSE){
    if (split){
        xx<-pivot_longer(x, c("precision", "recall"), names_to="type", values_to="kallisto")
        xx$coltojoin <- paste0(xx$tool, xx$splits, xx$type)
        yy<-pivot_longer(y, c("precision", "recall"), names_to="type", values_to="salmon")
        yy$coltojoin <- paste0(yy$tool, yy$splits, yy$type)
        

        aa <- left_join(xx,yy,by="coltojoin")
        if (!is.null(z)){
            zz<-pivot_longer(z, c("precision", "recall"), names_to="type", values_to="rsem")
            zz$coltojoin <- paste0(zz$tool, zz$splits, zz$type)
            aa <- left_join(aa,zz,by="coltojoin")
            g<-ggplot(aa)+
            geom_segment(aes(x=tool.x, xend=tool.x, y=kallisto, yend=salmon), color="grey")+
            geom_point(aes(x=tool.x, y=kallisto, color="#3300ff"), size=5)+
            geom_point(aes(x=tool.x, y=salmon, color="#f09c3b"), size=5)+
            geom_point(aes(x=tool.x, y=rsem, color="seagreen3"), size=5)+
            facet_grid(.~type.x)+coord_flip()
        } else {
            
            g <- ggplot(aa)+
                geom_segment(aes(x=tool.x, xend=tool.x, y=kallisto, yend=salmon), color="grey")+
                geom_point(aes(x=tool.x, y=kallisto), color="#3300ff", size=5)+
                geom_point(aes(x=tool.x, y=salmon), color="#f09c3b", size=5)+
                facet_grid(splits.x~type.x)+coord_flip()
        }
            
    } else {
        xx<-pivot_longer(x, c("precision", "recall"), names_to="type", values_to="kallisto")
        xx$coltojoin <- paste0(xx$tool, xx$type)
        yy<-pivot_longer(y, c("precision", "recall"), names_to="type", values_to="salmon")
        yy$coltojoin <- paste0(yy$tool, yy$type)
        
        aa <- left_join(xx,yy,by="coltojoin")
        if (!is.null(z)){
            zz<-pivot_longer(z, c("precision", "recall"), names_to="type", values_to="rsem")
            zz$coltojoin <- paste0(zz$tool, zz$type)
            
            aa <- left_join(aa,zz,by="coltojoin")
            g<-ggplot(aa)+
                geom_segment(aes(x=tool.x, xend=tool.x, y=kallisto, yend=salmon), color="grey")+
                geom_point(aes(x=tool.x, y=kallisto, color="Kallisto"), size=5)+
                geom_point(aes(x=tool.x, y=salmon, color="Salmon"), size=5)+
                geom_point(aes(x=tool.x, y=rsem, color="RSEM"), size=5)+
                facet_grid(.~type.x)+coord_flip()
        } else {
            g<-ggplot(aa)+
            geom_segment(aes(x=tool.x, xend=tool.x, y=kallisto, yend=salmon), color="grey")+
            geom_point(aes(x=tool.x, y=kallisto, color="Kallisto"), size=5)+
            geom_point(aes(x=tool.x, y=salmon, color="Salmon"), size=5)+
            facet_grid(.~type.x)+coord_flip()
        }
            
    }
    if (!is.null(zz)){
        g+theme_bw()+
            theme(axis.text=element_text(size=16), axis.title = element_text(size=20),
            strip.text.x = element_text(size = 20))+
            scale_color_manual(values=c("Kallisto"="#3300ff", "Salmon"="#f09c3b", "RSEM"="seagreen3"), 
                        name="Pseudoalignment tools",
                        breaks=c("Kallisto", "Salmon", "RSEM"))+
            xlab("Tools")
    } else {
        g+theme_bw()+
            theme(axis.text=element_text(size=16), axis.title = element_text(size=20),
            strip.text.x = element_text(size = 20))+
            scale_color_manual(values=c("Kallisto"="#3300ff", "Salmon"="#f09c3b"), 
                        name="Pseudoalignment tools",
                        breaks=c("Kallisto", "Salmon"))+
            xlab("Tools")
    }
            

}

compare_stager <- function(x){

}


### colors
pal_ramp <- function(values) {
  force(values)
  function(n) {
    if (n <= length(values)) {
      values[seq_len(n)]
    } else {
      colorRampPalette(values, alpha = TRUE)(n)
    }
  }
}

pal_adaptive <- function(name, palette, alpha = 1) {
  if (alpha > 1L | alpha <= 0L) stop("alpha must be in (0, 1]")

  raw_cols <- ggsci:::ggsci_db[[name]][[palette]]
  raw_cols_rgb <- col2rgb(raw_cols)
  alpha_cols <- rgb(
    raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
    alpha = alpha * 255L, names = names(raw_cols),
    maxColorValue = 255L
  )

  pal_ramp(unname(alpha_cols))
}

scale_color_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("colour", name, pal_adaptive(name, palette, alpha), ...)
}

scale_fill_adaptive <- function(name, palette, alpha = 1, ...) {
  ggplot2::discrete_scale("fill", name, pal_adaptive(name, palette, alpha), ...)
}