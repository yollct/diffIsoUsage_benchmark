thresholds <- data.frame(tools=c("iso_ktsp","drimseq","dexseq","dturtle","seqGSEA","cuffdiff","junctionseq","saturn","drimseq_stageR", "dexseq_stageR","saturn_stageR"),
                        thres=c(0.75,0.1,0.1,0.05,0.75,0.05,0.05,0.05,0.05,0.05,0.05))
row.names(thresholds) <- thresholds$tools

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
cal_pre_re <- function(methoddf, truthfile, thresholds, split=NULL, tx=TRUE){
    if (is.null(split)){
        mets <- colnames(methoddf)[2:ncol(methoddf)]
        xx<-1
        out <- do.call(cbind, sapply(mets, function(method){
            met <- methoddf[,method]
            names(met) <- methoddf$feature_id
            x<-thresholds[xx]
            #precisions <- sapply(thresholds, function(x){
            if (tx){
                tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
            } else {
                tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
            }
            pp <- length(names(met)[met<x])
            precisions <- tp/pp
            print(tp)
             #   })
            #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
            #recall <- sapply(thresholds, function(x){
            if (tx){
                tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
            } else {
                tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
            } 
            p <- length(truthfile[truthfile$status==1,]$feature_id)
            recall <- tp/p
            #    return(tp/p)})
            xx <- xx+1
            outout <- data.frame(fdr=c(precisions, recall), type=rep(c("precision", "recall"), each=length(precisions)), tool=rep(method, length(precisions)*2), thresholds=rep(x, 2))
            list(t(outout))
            }
        ))
    } else {
        out <- data.frame(row.names=c("fdr","type","tool","thresholds","splits"))
        for (i in unique(truthfile[split][,1])){
            splitted_truth <- truthfile %>% dplyr::filter(truthfile[split]==i)
            not_this_truth <- truthfile[truthfile[split]!=i,]$feature_id    
            splitted_method <- methoddf %>% dplyr::filter(!feature_id %in% not_this_truth)
            xx<-1
            thisout <- do.call(cbind, sapply(colnames(splitted_method)[2:ncol(splitted_method)], function(method){
                met <- splitted_method[,method]
                names(met) <- splitted_method$feature_id
                x <- thresholds[xx]
                #precisions <- sapply(thresholds, function(x){
                if (tx){
                    tp <- sum(names(met)[met<x] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                } else {
                    tp <- sum(names(met)[met<x] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                }
                pp <- length(names(met)[met<x])
                precisions <- tp/pp
                #})
                #print(truthfile$feature_id[!truthfile[truthfile$status==1,]$feature_id %in% names(met)[met<x]])
                #recall <- sapply(thresholds, function(x){
                if (tx){
                    tp <- sum(names(met)[met<x] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                    p <- length(splitted_truth[splitted_truth$status==1,]$feature_id)
                } else {
                    tp <- sum(names(met)[met<x] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                    p <- length(splitted_truth[splitted_truth$status==1,]$feature_id)
                }
            
                
                recall <- tp/p
            #})
                xx <- xx+1
                outout <- data.frame(fdr=c(precisions, recall), type=rep(c("precision", "recall"), each=length(precisions)), tool=rep(method, length(precisions)*2), thresholds=rep(x,2), splits=rep(i, length(precisions)*2))
                list(t(outout))
            }))
            out <- cbind(out, thisout)
        }
    }

    return(out)
}

cal_pre_re_curve <- function(methoddf, truthfile, split=NULL){
    if (is.null(split)){
        mets <- colnames(methoddf)[2:ncol(methoddf)]
        out <- do.call(rbind, sapply(mets, function(method){
            met <- methoddf[,method]
            names(met) <- methoddf$feature_id
            
            thresholds <- seq(0,max(met), 0.1)
            precisions <- sapply(thresholds, function(x){
                tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
                pp <- length(names(met)[met<x])
                tp/pp
                })

                
            recall <- sapply(thresholds, function(x){
                tp <- sum(names(met)[met<x] %in% truthfile[truthfile$status==1,]$feature_id) 
                p <- length(truthfile[truthfile$status==1,]$feature_id)
                return(tp/p)})

            outout <- data.frame(fdr=c(precisions, recall), type=rep(c("precision", "recall"), each=length(precisions)), tool=rep(method, length(precisions)*2), thresholds=rep(thresholds, 2))
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
                met <- splitted_method[,method]
                names(met) <- splitted_method$feature_id
                precisions <- sapply(thresholds, function(x){
                    tp <- sum(names(met)[met<x] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                    pp <- length(names(met)[met<x])
                    tp/pp
                })
                
                recall <- sapply(thresholds, function(x){
                    tp <- sum(names(met)[met<x] %in% splitted_truth[splitted_truth$status==1,]$feature_id)
                    p <- length(splitted_truth[splitted_truth$status==1,]$feature_id)
                    tp/p
                })
                outout <- data.frame(fdr=c(precisions, recall), type=rep(c("precision", "recall"), each=length(precisions)), tool=rep(method, length(precisions)*2), thresholds=rep(thresholds,2), splits=rep(i, length(precisions)*2))
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
        wideoutputpr$splits <- factor(wideoutputpr$splits, level=c("DTE","DTU","DTEDTU"))
        wideoutputpr$recall[is.nan(wideoutputpr$recall)] = 0
    }

    return(wideoutputpr)
}

compare_tools <- function(x,y,split=FALSE){
    if (split){
        xx<-pivot_longer(x, c("precision", "recall"), names_to="type", values_to="kallisto")
        xx$coltojoin <- paste0(xx$tool, xx$splits, xx$type)
        yy<-pivot_longer(y, c("precision", "recall"), names_to="type", values_to="salmon")
        yy$coltojoin <- paste0(yy$tool, yy$splits, yy$type)

        aa <- left_join(xx,yy,by="coltojoin")
        g <- ggplot(aa)+
            geom_segment(aes(x=tool.x, xend=tool.x, y=kallisto, yend=salmon), color="grey")+
            geom_point(aes(x=tool.x, y=kallisto), color="#3300ff", size=5)+
            geom_point(aes(x=tool.x, y=salmon), color="#f09c3b", size=5)+
            facet_grid(splits.x~type.x)+coord_flip()
            
    } else {
        xx<-pivot_longer(x, c("precision", "recall"), names_to="type", values_to="kallisto")
        xx$coltojoin <- paste0(xx$tool, xx$type)
        yy<-pivot_longer(y, c("precision", "recall"), names_to="type", values_to="salmon")
        yy$coltojoin <- paste0(yy$tool, yy$type)

        aa <- left_join(xx,yy,by="coltojoin")
        g<-ggplot(aa)+
            geom_segment(aes(x=tool.x, xend=tool.x, y=kallisto, yend=salmon), color="grey")+
            geom_point(aes(x=tool.x, y=kallisto, color="Kallisto"), size=5)+
            geom_point(aes(x=tool.x, y=salmon, color="Salmon"), size=5)+
            facet_grid(.~type.x)+coord_flip()
            
    }
    g+theme_bw()+
        theme(axis.text=element_text(size=16), axis.title = element_text(size=20),
        strip.text.x = element_text(size = 20))+
        scale_color_manual(values=c("Kallisto"="#3300ff", "Salmon"="#f09c3b"), 
                    name="Pseudoalignment tools",
                    breaks=c("Kallisto", "Salmon"))+
        xlab("Tools")
        

}