library(stringr)
library(grDevices)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(dplyr)
library(reshape)
library(ggplot2)
library(mixdist)
#######################################################################
########################  STEP1  ######################################

reformat <- function(MID_unique,Rhapsody_input,peptidelist){
  miss <- setdiff(peptidelist$name,colnames(MID_unique))
  MID_unique[miss] <- 0
  MID_unique <- MID_unique[peptidelist$name]
  MID_unique <- MID_unique[intersect(rownames(Rhapsody_input),rownames(MID_unique)),]
  miss <- setdiff(rownames(Rhapsody_input),rownames(MID_unique))
  MID_unique[miss,] <- 0
  MID_unique <- MID_unique[rownames(Rhapsody_input),]
}

#######################################################################
########################  STEP2  ######################################

#functions
specificity <- function(counts){
  if(max(counts)==0){
    return("NEG")
  }else{
    counts <- sort(counts,decreasing = T)
    rawdiff <- diff(counts)/diff(seq(length(counts)))
    inflection <- max(which(rawdiff == min(rawdiff, na.rm=TRUE)))
    spec <- names(counts[1:inflection])
    return(paste(spec,collapse = "|"))
  }
}

specificity_MIDs <- function(counts){
  counts <- sort(counts,decreasing = T)
  rawdiff <- diff(counts[1:8])/diff(seq(8))
  inflection <- max(which(rawdiff == min(rawdiff[1:8], na.rm=TRUE)))
  MIDs <- counts[1:inflection]
  return(paste(MIDs,collapse = "|"))
}


##Main function
tetramer <- function(MID_unique, peptidelist, pep_dist_thr = 3, breaks = 50, asinh_scale=4, asinh_scale2=4){
  cat = levels(as.factor(peptidelist$category))[!levels(as.factor(peptidelist$category))=="Empty"]
  summary <- aggregate(. ~ category, data=data.frame(category=peptidelist$category,t(MID_unique),check.names = F),FUN = sum)
  rownames(summary) <- summary$category
  summary <- as.data.frame(t(summary[,-1]))
  
  if(length(cat) == 1){
    pdf(file="figures/QC_negative_threshold.pdf",width = 4,height=4)
    cat.name <- peptidelist[peptidelist$category == cat,]$name
    hist <- hist(asinh(summary[,cat]/asinh_scale),breaks = breaks)
    df <- data.frame(mid=hist$mids,cou=hist$counts)
    for(i in seq(0.1,5,by=0.1)){
      res <- try(mix.model <- mix(as.mixdata(df),mixparam(mu=c(0,2),sigma=i)))
      if(inherits(res, "try-error")){
        next
      }else{
        break
      }
    }
    f <- function(x){dnorm(x,res$parameters$mu[1],res$parameters$sigma[1]) * res$parameters$pi[1] - dnorm(x,res$parameters$mu[2],res$parameters$sigma[2]) * res$parameters$pi[2]}
    thr <- uniroot(f,interval = c(res$parameters$mu[1],res$parameters$mu[2]))$root
    cutoff <- ceiling((0.5*exp(-thr)*(exp(2*thr)-1))*asinh_scale)

    #cat.d <- density(log10(summary[,cat]+1),bw = 0.1)
    #local.min.cat <- optimize(approxfun(cat.d$x,cat.d$y),interval=c(0,1.5))$minimum
    #cutoff.cat <- round(10^(local.min.cat))
    #hist(log10(summary[,cat]+1),breaks = 50,probability = T,col = "green",main = cat,xlab="log10MID",border = "green")
    plot(res,main = cat,xlab="asinh of MID")
    #lines(cat.d, col="red",lwd=2,lty=2)
    abline(v=thr, col=c("blue"), lty=2,lwd=2)
    neg <- cutoff
    names(neg) <- cat
    dev.off()
    
  }else{
    cat1 = as.character(cat[1])
    cat2 = as.character(cat[2])
    pdf(file="figures/QC_negative_threshold.pdf",width = 8,height=8)
    layout(matrix(c(1,2,3,4),2,2,byrow = T))
    ##summary QC density scatter plot:
    x <- densCols(log10(summary[,cat1]+1),log10(summary[,cat2]+1),colramp=colorRampPalette(c("black", "white")))
    plot(log10(summary[,cat1]+1),log10(summary[,cat2]+1),pch=20,col=colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)[col2rgb(x)[1,] + 1L],xlab=paste("log10MID count of", cat1, "peptides"),ylab=paste("log10MID count of", cat2, "peptides"))
    plot(log10(summary[,cat1]+1),log10(summary[,cat2]+1),pch=20,col=colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)[col2rgb(x)[1,] + 1L],xlab=paste("log10MID count of", cat1, "peptides"),ylab=paste("log10MID count of", cat2, "peptides"),xlim=c(0,2.5),ylim=c(0,2.5))
    cat1.name <- peptidelist[peptidelist$category == cat1,]$name
    cat2.name <- peptidelist[peptidelist$category == cat2,]$name
  
    ##identify negative threshold:
    his <- (summary[,cat1] < 3 & summary[,cat2] > 0) | (summary[,cat2] < 3 & summary[,cat1] > 0)
    ###cat1
    hist.cat1 <- hist(asinh(summary[his, cat1]/asinh_scale),breaks = breaks,plot=FALSE)
    #hist.cat1 <- hist(asinh(summary[,cat1]/asinh_scale),breaks = breaks,plot=FALSE)
    df.cat1 <- data.frame(mid=hist.cat1$mids,cou=hist.cat1$counts)
    for(i in seq(0.1,5,by=0.1)){
      print(i)
      res.cat1 <- try(mix.model <- mix(as.mixdata(df.cat1),mixparam(mu=c(0,2),sigma=i)))
      if(inherits(res.cat1, "try-error")){
        next
      }else{
        break
      }
    }
    f.cat1 <- function(x){dnorm(x,res.cat1$parameters$mu[1],res.cat1$parameters$sigma[1]) * res.cat1$parameters$pi[1] - dnorm(x,res.cat1$parameters$mu[2],res.cat1$parameters$sigma[2]) * res.cat1$parameters$pi[2]}
    thr.cat1 <- uniroot(f.cat1,interval = c(res.cat1$parameters$mu[1],res.cat1$parameters$mu[2]))$root
    cutoff.cat1 <- ceiling((0.5*exp(-thr.cat1)*(exp(2*thr.cat1)-1))*asinh_scale)
#    cat1.d <- density(log10(summary[summary[,cat1] >= summary[,cat2],cat1]+1),bw = density.bw)
#    local.min.cat1 <- optimize(approxfun(cat1.d$x,cat1.d$y),interval=c(0,1.5))$minimum
#    cutoff.cat1 <- round(10^(local.min.cat1))

    ###cat2
    hist.cat2 <- hist(asinh(summary[his, cat2]/asinh_scale2),breaks = breaks,plot=FALSE)
    #hist.cat2 <- hist(asinh(summary[,cat2]/asinh_scale2),breaks = breaks,plot=FALSE)
    df.cat2 <- data.frame(mid=hist.cat2$mids,cou=hist.cat2$counts)
    for(i in seq(0.1,5,by=0.1)){
      print(i)
      res.cat2 <- try(mix.model <- mix(as.mixdata(df.cat2),mixparam(mu=c(0,2),sigma=i)))
      if(inherits(res.cat2, "try-error")){
        next
      }else{
        break
      }
    }
    f.cat2 <- function(x){dnorm(x,res.cat2$parameters$mu[1],res.cat2$parameters$sigma[1]) * res.cat2$parameters$pi[1] - dnorm(x,res.cat2$parameters$mu[2],res.cat2$parameters$sigma[2]) * res.cat2$parameters$pi[2]}
    thr.cat2 <- uniroot(f.cat2,interval = c(res.cat2$parameters$mu[1],res.cat2$parameters$mu[2]))$root
    cutoff.cat2 <- ceiling((0.5*exp(-thr.cat2)*(exp(2*thr.cat2)-1))*asinh_scale2)
#    cat2.d <- density(log10(summary[summary[,cat1] <= summary[,cat2],cat2]+1),bw = density.bw)
#    local.min.cat2 <- optimize(approxfun(cat2.d$x,cat2.d$y),interval=c(0,1.5))$minimum
#    cutoff.cat2 <- round(10^(local.min.cat2))

    neg <- c(cutoff.cat1,cutoff.cat2)
    names(neg) <- c(cat1,cat2)
    
#    hist(log10(summary[summary[,cat1] > summary[,cat2],cat1]+1),breaks = 50,probability = T,col = "green",main = cat[1],xlab="log10MID",border = "green")
#    lines(cat1.d, col="red",lwd=2,lty=2)
    plot(res.cat1,main = cat[1],xlab="asinh of MID")
    abline(v=thr.cat1, col=c("blue"), lty=2,lwd=2)
    text(thr.cat1,0.2,print(paste0("MID cutoff =", cutoff.cat1)),col="red")
#    hist(log10(summary[summary[,cat1] < summary[,cat2],cat2]+1),breaks = 50,probability = T,col = "green",main = cat[2],xlab="log10MID",border = "green")
#    lines(cat2.d, col="red",lwd=2,lty=2)
    plot(res.cat2,main = cat[2],xlab="asinh of MID")    
    abline(v=thr.cat2, col=c("blue"), lty=2,lwd=2)
    text(thr.cat2,0.2,print(paste0("MID cutoff =", cutoff.cat2)),col="red")
    dev.off()
  
  }
  
  spec.summary <- as.data.frame(matrix(nrow = nrow(MID_unique),ncol = 7))
  colnames(spec.summary) <- c("totalMID","Top1MID","Top1peptide","specificity","specificity.MIDs","SignalRatio","gene")
  rownames(spec.summary) <- rownames(MID_unique)
  spec.summary$totalMID <- rowSums(MID_unique)
  spec.summary$Top1MID <- apply(MID_unique,1,max)
  spec.summary$Top1peptide <- ifelse(spec.summary$Top1MID == 0, "NEG", colnames(MID_unique)[apply(MID_unique,1,which.max)])
  inflection <- apply(as.matrix(MID_unique),1,function(x) specificity(x))
  spec.summary$specificity.MIDs <- apply(as.matrix(MID_unique),1,function(x) specificity_MIDs(x))
  
  seq <- peptidelist$AAsequence
  names(seq) <- peptidelist$name
  
  if(length(cat) == 1){
    for(i in 1:nrow(spec.summary)){
      count.rank <- order(MID_unique[i,],decreasing = T)
      lv <- stringdist::stringdist(seq[colnames(MID_unique)[count.rank[1]]],seq[colnames(MID_unique)[count.rank[2:length(count.rank)]]],method = "lv")
      if(spec.summary[i,"Top1MID"] < neg+1){
        spec.summary[i,"specificity"] <- "NEG"
        spec.summary[i,"specificity.MIDs"] <- spec.summary[i,"Top1MID"]
      }else if(lv[1] <= pep_dist_thr){
        cross.num <- min(which(lv > pep_dist_thr))
        spec.summary[i,"specificity"] <- paste(colnames(MID_unique)[count.rank[1:cross.num]],collapse ="|")
        spec.summary[i,"specificity.MIDs"] <- paste(MID_unique[i, count.rank[1:cross.num]],collapse ="|")
      }else if(grepl("\\|",inflection[i])){
        spec.summary[i,"specificity"] <- paste(inflection[i],"undetermined", sep = ";")
      }else{
        spec.summary[i,"specificity"] <- inflection[i]
      }
    }
  }else{
    for(i in 1:nrow(spec.summary)){
      count.rank <- order(MID_unique[i,],decreasing = T)
      lv <- stringdist::stringdist(seq[colnames(MID_unique)[count.rank[1]]],seq[colnames(MID_unique)[count.rank[2:length(count.rank)]]],method = "lv")
      if(spec.summary[i,"Top1MID"] < neg[1]+1 & spec.summary[i,"Top1peptide"] %in% c(cat1.name,"NEG")){
        spec.summary[i,"specificity"] <- "NEG"
        spec.summary[i,"specificity.MIDs"] <- spec.summary[i,"Top1MID"]
      }else if(spec.summary[i,"Top1MID"] < neg[2]+1 & spec.summary[i,"Top1peptide"] %in% c(cat2.name,"NEG")){
        spec.summary[i,"specificity"] <- "NEG"
        spec.summary[i,"specificity.MIDs"] <- spec.summary[i,"Top1MID"]
      }else if(lv[1] <= pep_dist_thr){
        cross.num <- min(which(lv > pep_dist_thr))
        spec.summary[i,"specificity"] <- paste(colnames(MID_unique)[count.rank[1:cross.num]],collapse ="|")
        spec.summary[i,"specificity.MIDs"] <- paste(MID_unique[i, count.rank[1:cross.num]],collapse ="|")
      }else if(grepl("\\|",inflection[i])){
        spec.summary[i,"specificity"] <- paste(inflection[i],"undetermined", sep = ";")
      }else{
        spec.summary[i,"specificity"] <- inflection[i]
      }
    }
  }
  
  gene <- peptidelist$gene
  names(gene) <- peptidelist$name
  spec.summary$gene <- gene[spec.summary$specificity]
  cr <- !grepl("undetermined",spec.summary$specificity) & grepl("\\|",spec.summary$specificity)
  spec.summary[cr,"gene"] <- "crossreactive"
  signal <- sapply(spec.summary$specificity.MIDs, function(x){
    return(sum(as.numeric(unlist(strsplit(x,split = "\\|")))))
  })
  spec.summary$SignalRatio <- signal / spec.summary$totalMID
  
  sorted.matrix <- t(apply(MID_unique,1,sort,decreasing=T))  
  if(length(cat)==1){
    sorted.cat <- subset(sorted.matrix, spec.summary$specificity %in% cat.name)
    sorted.cat.melt <- melt(sorted.cat[head(order(sorted.cat[,1],decreasing = T),100),1:8])
    p <- ggplot(sorted.cat.melt,aes(X2,value)) + geom_line(aes(group=X1),color="deepskyblue") + theme_bw() +labs(x="peptide rank",y="MID counts",title="mono-specific peptide MID ranks") + scale_x_continuous(breaks=c(1,3,5,7,9))
    
    if(sum(cr) > 0){
      sorted.cross <- subset(sorted.matrix,cr)
      if(sum(cr) > 100){
        sorted.cross.melt <- melt(sorted.cross[head(order(sorted.cross[,1],decreasing = T),100),1:8])
      }else{
        sorted.cross.melt <- melt(sorted.cross[,1:8])
      }
      p1 <- ggplot(sorted.cross.melt,aes(X2,value)) + geom_line(aes(group=X1),color="deepskyblue") + theme_bw() +labs(x="peptide rank",y="MID counts",title="cross-reactive peptide MID ranks") + scale_x_continuous(breaks=c(1,3,5,7,9))
    }
      
    if(sum(cr) > 0){
      ggsave(file="figures/QC_peprank.pdf", gridExtra::arrangeGrob(p,p1,nrow = 2), width = 4, height=8)
    }else{
      ggsave(file="figures/QC_peprank.pdf", p, width = 4, height=4)
    }
  }else{
    sorted.cat1 <- subset(sorted.matrix, spec.summary$specificity %in% cat1.name)
    sorted.cat1.melt <- melt(sorted.cat1[head(order(sorted.cat1[,1],decreasing = T),100),1:8])
    p1 <- ggplot(sorted.cat1.melt,aes(X2,value)) + geom_line(aes(group=X1),color="deepskyblue") + theme_bw() +labs(x="peptide rank",y="MID counts",title=paste(cat1,"peptide MID ranks")) + scale_x_continuous(breaks=c(1,3,5,7,9))
    
    sorted.cat2 <- subset(sorted.matrix, spec.summary$specificity %in% cat2.name)
    sorted.cat2.melt <- melt(sorted.cat2[head(order(sorted.cat2[,1],decreasing = T),100),1:8])
    p2 <- ggplot(sorted.cat2.melt,aes(X2,value)) + geom_line(aes(group=X1),color="deepskyblue") + theme_bw() +labs(x="peptide rank",y="MID counts",title=paste(cat2,"peptide MID ranks")) + scale_x_continuous(breaks=c(1,3,5,7,9))
    
    if(sum(cr) > 0){
      sorted.cross <- subset(sorted.matrix,cr)
      if(sum(cr) > 100){
        sorted.cross.melt <- melt(sorted.cross[head(order(sorted.cross[,1],decreasing = T),100),1:8])
      }else{
        sorted.cross.melt <- melt(sorted.cross[,1:8])
      }
      p3 <- ggplot(sorted.cross.melt,aes(X2,value)) + geom_line(aes(group=X1),color="deepskyblue") + theme_bw() +labs(x="peptide rank",y="MID counts",title="cross-reactive peptide MID ranks") + scale_x_continuous(breaks=c(1,3,5,7,9))
    }
    
    if(sum(cr) > 0){
      ggsave(file="figures/QC_peprank.pdf", gridExtra::arrangeGrob(p1,p2,p3,nrow = 2))
    }else{
      ggsave(file="figures/QC_peprank.pdf",gridExtra::arrangeGrob(p1,p2,nrow = 2))
    }
  }
  spec.summary[is.na(spec.summary)] <- ""
  out <- list("spectable"=spec.summary,"MIDsummary"=summary,"negative_threshold"=neg)
  return(out)
}


#######################################################################
########################  STEP3  ######################################

###custom colorramps function:
# Returns a vector of 'num.colors.in.palette'+1 colors. The first 'cutoff.fraction'
# fraction of the palette interpolates between colors[1] and colors[2], the remainder
# between colors[3] and colors[4]. 'num.colors.in.palette' must be sufficiently large
# to get smooth color gradients.
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 3)
  ramp1 <- colorRampPalette(colors[1:1])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[2:3])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}

memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}

tcrcorrect <- function(tetTCRtable){
  tcr <- tetTCRtable$TCR
  tcr.spec.correct <- sapply(tcr, function(x) {
    tp <- table(tetTCRtable[tetTCRtable$TCR == x,]$specificity)
    if(length(tp)==1){
      return(names(tp))
    }else{
      tp1 <- tp[!grepl("undetermined|NEG",names(tp))]
      if(length(tp1)==0){
        return(names(which.max(tp)))
      }else if(length(tp1)==1){
        if(tp1 > 1){
          return(names(tp1))
        }else{
          return(names(which.max(tp)))
        }
      }else{
        top1 <- sort(tp1,decreasing = T)[1]
        top2 <- sort(tp1,decreasing = T)[2]
        if(top1 > 2*top2){
          return(names(top1))
        }else{
          return("undetermined")
        }
      }
    }
  })
  return(tcr.spec.correct)
}

###tetramer+TCR analysis:
tetTCR <- function(spectable,TCRtable,peptidelist){
  TCRtable$TCR <- paste(TCRtable$TCRa_1st_cdr_aa,TCRtable$TCRb_1st_cdr_aa,sep="_")
  tmp <- cbind(spectable,TCRtable)
  clean <- !grepl("undetermined|NEG|HCV",tmp$specificity) & tmp[,"TCRa_1st_cdr_aa"] != "" & tmp[,"TCRb_1st_cdr_aa"] != ""
  crosstab <- table(tmp[clean,]$TCR,tmp[clean,]$specificity)
  top10.crosstab <- crosstab[head(order(rowSums(crosstab),decreasing = T),10),]
  top10.crosstab <- top10.crosstab[,colSums(top10.crosstab)>0]
  cols <- makeColorRampPalette(c("white", "white", "red"),1/ max(crosstab),100)
  pheatmap(t(log10(memoSort(top10.crosstab)+1)),cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,color = cols,border_color = NA,
           filename = "figures/Top10TCR_specificities.pdf",width = 6,height=3, main="top10TCR X Specificity")
  
#  uniq.pair.idx <- rowSums(tmp[,c("TCRb_1st_cdr_aa","TCRb_2nd_cdr_aa")] != "")==1 & rowSums(tmp[,c("TCRa_1st_cdr_aa","TCRa_2nd_cdr_aa")] != "")==1 & tmp[,"gene"] != "crossreactive"
#  multi.idx <- rowSums(tmp[,c("TCRb_1st_cdr_aa","TCRb_2nd_cdr_aa")] != "")==2 & grepl("undetermined",tmp$specificity)
#  rest.idx <- !(uniq.pair.idx | multi.idx)
#  correct <- tcrcorrect(tmp[uniq.pair.idx,])
#   
#  final <- spectable
#  final[c("correct.status","specificity.final","gene.final")] <- ""
#  final[uniq.pair.idx,"specificity.final"] <- correct
#  final[uniq.pair.idx,"correct.status"] <- "correctable"
#  final[multi.idx,"specificity.final"] <- "multiplets"
#  final[multi.idx,"correct.status"] <- "multiplets"
#  final[rest.idx,"specificity.final"] <- final[rest.idx,]$specificity
#  final[rest.idx,"correct.status"] <- "uncorrectable"
#  
#  gene <- peptidelist$gene
#  names(gene) <- peptidelist$name
#  final$gene.final <- gene[final$specificity.final]
#  cr <- !grepl("undetermined",final$specificity.final) & grepl("\\|",final$specificity.final)
#  final[cr,"gene.final"] <- "crossreactive"
#  
#  complete <- cbind(final, TCRtable)
#  complete$specificity.frequency <- sapply(complete$specificity.final, function(x) sum(complete$specificity.final == x))
#  complete$tcr.frequency <- sapply(complete$TCR, function(x) sum(complete$TCR == x))
#
#  summary <- data.frame(row.names = c("number of cells","number of cells with single TCR chains","number of cells with TCRa","number of cells with TCRb",
#                                      "number of cells with unique specificity","number of multiplets","specificity fdr"))
#  total <- nrow(complete)
#  uniqTCR <- sum(rowSums(tmp[,c("TCRb_1st_cdr_aa","TCRb_2nd_cdr_aa")] != "")==1 & rowSums(tmp[,c("TCRa_1st_cdr_aa","TCRa_2nd_cdr_aa")] != "")==1)
#  a <- sum(complete$TCRa_1st_cdr_aa != "")
#  b <- sum(complete$TCRb_1st_cdr_aa != "")
#  uniqspec <- sum(!grepl("undetermined|NEG|multiplets",final$specificity.final))
#  multiplets <- sum(multi.idx)
#  corrected <- final[final$correct.status=="correctable",]
#  ts <- !grepl("undetermined|NEG",corrected$specificity.final)
#  diff <- sum(corrected[ts,"specificity.final"] != corrected[ts,"specificity"])
#  fdr <- diff/sum(ts)
#  summary$stat <- c(as.character(c(total,uniqTCR,a,b,uniqspec,multiplets)),fdr)
  
#  out <- list("final.specificity"=final,"completetable"=complete,"summary"=summary)
#  return(out)
}


#######################################################################
########################  STEP5  ######################################

### Analysis of differences of specificities

spec_analysis <- function(tSNE_out,Rhapsody_input,peptidelist){
  spec <- unique(tSNE_out$spec.corrected)[!unique(tSNE_out$spec.corrected) %in% c("NEG","undetermined")]
  spec.frequency <- sapply(spec, function(x) sum(combo0920_tsne$spec.corrected==x))
  spec <- spec[spec.frequency>5]
  pdf(file="figures/specificity_PCAplots.pdf",width = 6,height=6)
  
  ###### V1 ######
  ###### calculate the percentage of specificity in each cluster ---overall of all donors#####
  tSNE_out.filter <- tSNE_out[tSNE_out$spec.corrected %in% spec,]
  spec.pheno.count <- table(tSNE_out.filter$spec.corrected,tSNE_out.filter$cluster)
  spec.pheno.percent <- spec.pheno.count/rowSums(spec.pheno.count)
  spec.pheno.pca <- prcomp(spec.pheno.percent)
  spec.pheno.out <- as.data.frame(spec.pheno.pca$x[,1:2])
  pc1.var <- summary(spec.pheno.pca)$importance[2,1] * 100
  pc2.var <- summary(spec.pheno.pca)$importance[2,2] * 100
  spec.pheno.out$category <- cross_assign_category(rownames(spec.pheno.out),peptidelist)
  spec.pheno.out$gene <- cross_assign_gene(rownames(spec.pheno.out),peptidelist)
  spec.pheno.out$frequency <- sapply(rownames(spec.pheno.out), function(x) sum(tSNE_out$spec.corrected==x))
  print(ggplot(spec.pheno.out,aes(PC1,PC2))+geom_point(aes(color=category))
        +labs(x=paste("PC1 (",pc1.var,"%)",sep = ""),y=paste("PC2 (",pc2.var,"%)",sep = ""),title="PCA_v1")+theme_bw()) 
  
  ###### V2 ######
  ###### calculate the percentage of specificity in each cluster ---breakdown to donors#####
  cell.remove <- grepl("negative",tSNE_out$Sample_Tag) | tSNE_out$Sample_Tag %in% c("Multiplet","Undetermined") | !(tSNE_out$spec.corrected %in% spec)
  spec.donor <- tSNE_out[!cell.remove,] %>% count(spec.corrected,Sample_Tag,cluster)
  spec.donor.percent <- spec.donor %>% group_by(spec.corrected,Sample_Tag) %>% mutate(prop= n / sum(n))
  spec.donor.spread <- spec.donor.percent %>% spread(cluster, prop)
  spec.donor.spread[is.na(spec.donor.spread)] <- 0
  spec.donor.pca <- prcomp(spec.donor.spread[,3:ncol(spec.donor.spread)])
  spec.donor.pca.out <- as.data.frame(spec.donor.pca$x[,1:2])
  spec.donor.pca.out$donor <- spec.donor.spread$Sample_Tag
  spec.donor.pca.out$specificity <- spec.donor.spread$spec.corrected
  spec.donor.pca.out$category <- cross_assign_category(spec.donor.pca.out$specificity,peptidelist)
  spec.donor.pca.out$gene <- cross_assign_gene(spec.donor.pca.out$specificity,peptidelist)
  pc1.var <- summary(spec.donor.pca)$importance[2,1] * 100
  pc2.var <- summary(spec.donor.pca)$importance[2,2] * 100
  print(ggplot(spec.donor.pca.out,aes(PC1,PC2))+geom_point(aes(color=category))
        +labs(x=paste("PC1 (",pc1.var,"%)",sep = ""),y=paste("PC2 (",pc2.var,"%)",sep = ""),title="PCA_v2")+theme_bw())
  print(ggplot(spec.donor.pca.out, aes(PC1, PC2, colour = spec.corrected)) + scale_color_brewer(palette = "Paired") + 
          geom_point(data = spec.donor.pca.out[,1:2], colour = "grey", alpha = .2) + 
          facet_wrap(~ gene) + guides(colour = FALSE) + theme_bw())
  
  ###### V3 ######
  ###### calculate 90% percentile expression of each gene in specificity*cluster combination
  cluster <- levels(tSNE_out$cluster)
  gene <- colnames(Rhapsody_input)
  quantile_expr <- sapply(spec, function(x){
    sapply(gene, function(y){
      sapply(cluster, function(z){
        index <- which(tSNE_out$spec.corrected == x & tSNE_out$cluster == z)
        expr <- quantile(Rhapsody_input[index,y],0.9)
        return(expr)
      })
    })
  })
  colnames(quantile_expr) <- spec
  rownames(quantile_expr) <- apply(format(expand.grid(cluster,gene)),1,paste,collapse="_")
  quantile_expr[is.na(quantile_expr)] <- 0
  quantile_expr <- quantile_expr[rowSums(quantile_expr) > 0,]
  pca_v3 <- prcomp(t(log2(quantile_expr+1)))
  pca_v3.out <- as.data.frame(pca_v3$x[,1:2])
  pc1.var <- summary(pca_v3)$importance[2,1] * 100
  pc2.var <- summary(pca_v3)$importance[2,2] * 100
  pca_v3.out$category <- cross_assign_category(rownames(pca_v3.out),peptidelist)
  pca_v3.out$gene <- cross_assign_gene(rownames(pca_v3.out),peptidelist)
  pca_v3.out$frequency <- sapply(rownames(pca_v3.out), function(x) sum(tSNE_out$spec.corrected==x))
  print(ggplot(pca_v3.out,aes(PC1,PC2))+geom_point(aes(color=category)) 
        +labs(x=paste("PC1 (",pc1.var,"%)",sep = ""),y=paste("PC2 (",pc2.var,"%)",sep = ""),title="PCA_v3")+theme_bw())
  dev.off()
  return(list(specificity_summary = spec.donor.percent, pca_v1 = spec.pheno.out, pca_v2 = spec.donor.pca.out,
              quantile_expr = quantile_expr,pca_v3 = pca_v3.out))
}


