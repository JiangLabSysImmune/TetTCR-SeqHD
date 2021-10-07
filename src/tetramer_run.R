source('./src/tetramer_rev4.R')
library(gridExtra)
library(colorRamps)
####run commands####
args <- commandArgs(TRUE)
rhapsody <- read.csv(args[1], row.names = 1, comment.char = "#",check.names=F)
tetramer_dbec <- read.csv("./out/tetramer_dbec.csv", row.names=1,check.names = F)
peptide <- read.csv("./out/peptide.csv", stringsAsFactors = F)
tcra <- read.csv("./out/tcra_final.csv",fill = TRUE,header = T,row.names = 1,sep=";",check.names=F)
tcra[is.na(tcra)] <- ""
tcrb <- read.csv("./out/tcrb_final.csv",fill = TRUE,header = T,row.names = 1,sep=";",check.names=F)
tcrb[is.na(tcrb)] <- ""
sampletag <- read.csv("./out/sampletag.csv",header = T, row.names = 1)
cluster <- read.csv("./out/obs.csv",row.names = 1)
umapcoord <- read.csv("./out/obsm.csv",row.names = rownames(cluster))
#tsnecoord <- read.csv(args[2], row.names = rownames(rhapsody), comment.char = "#")
ggplot(umapcoord, aes(X_umap1,X_umap2)) + geom_point(aes(color=as.factor(cluster$leiden))) + labs(colour = "Cluster") + theme_bw()
ggsave("./figures/umap_cluster.pdf",width=5,height=4)
ggplot(umapcoord, aes(X_umap1,X_umap2)) + geom_point(aes(color=as.factor(sampletag$tag.final))) + scale_color_brewer(palette = "Set1") + labs(colour = "SampleTag") + theme_bw()
ggsave("./figures/umap_sampletag.pdf",width=6,height=4)

tcr <- cbind(tcra,tcrb)
colnames(tcr) <- c(paste("TCRa",colnames(tcra),sep = "_"), paste("TCRb",colnames(tcrb),sep = "_"))
rownames(tcr) <- rownames(rhapsody)
tetramer.reformat <- reformat(tetramer_dbec,rhapsody,peptide)
tetramer.out <- tetramer(tetramer.reformat, peptide,breaks=50,asinh_scale=4,asinh_scale2=4)
tetTCR(tetramer.out$spectable,tcr,peptide)
ab.gene <- colnames(rhapsody)[grep("pAbO", colnames(rhapsody))]
rhapsody.scale <- scale(log10(rhapsody+1))
p <- lapply(ab.gene, function(x) 
	ggplot(umapcoord, aes(X_umap1,X_umap2))
	+ geom_point(aes(color=as.numeric(rhapsody.scale[,x])), show.legend = F, size = 0.1)
	+ scale_color_gradientn(colours = matlab.like2(100))
	+ labs(title=x) + theme_bw() + theme(plot.title = element_text(size=10)))
ggsave("./figures/AbSeq_expression.pdf",marrangeGrob(p, nrow=3, ncol=3))

#tsne.plots <- tsnecoord
#tsne.plots$spec.final <- specificity.correct$completetable$specificity.final
#tsne.plots$sampletag <- sampletag$tag.final
#tsne.plots$cluster <- cluster$leiden
#p1 <- ggplot(tsne.plots[specificity.correct$completetable$specificity.frequency >= 10,], aes(x = Coordinate_1, y = Coordinate_2, colour = spec.final)) + geom_point(data = tsne.plots[,1:2], colour = "grey", alpha = .2, size = 0.1) + labs(title = "specificities on tSNE|>=10 cells|donors together") + geom_point(size = 0.3) + facet_wrap(~ spec.final) + guides(colour = FALSE) + theme_bw() + theme(strip.text = element_text(size = 5))
#p2 <- ggplot(tsne.plots[specificity.correct$completetable$specificity.frequency >= 10,], aes(x = Coordinate_1, y = Coordinate_2, colour = sampletag)) + scale_color_brewer(palette = "Set1") + geom_point(data = tsne.plots[,1:2], colour = "grey", alpha = .2, size = 0.1) + labs(title = "specificities on tSNE|>=10 cells|donors seperate") + guides(colour = guide_legend(override.aes = list(size=5))) + geom_point(size = 0.3) + facet_wrap(~ spec.final) + theme_bw() + theme(strip.text = element_text(size = 5))
#ggsave("./figures/specificity_tsne_plots.png", marrangeGrob(list(p1,p2),nrow=1,ncol=1),device="png",dpi=300)
#spec.filter <- specificity.correct$completetable$specificity.frequency >= 10 & !grepl("undetermined|NEG",tsne.plots$spec.final)
#cont.table <- table(tsne.plots[spec.filter,]$spec.final, tsne.plots[spec.filter,]$cluster)
#cont.table <- cont.table/rowSums(cont.table)
#pheatmap(cont.table,color=colorRampPalette(c("purple","yellow"))(100),border_color = NA, filename = "figures/specificity x cluster.pdf")

#tsne.plots$TCR <- specificity.correct$completetable$TCR
#tmp <- specificity.correct$completetable
#tmp.unique.tcr <- rowSums(tmp[,c("TCRb_1st_cdr_aa","TCRb_2nd_cdr_aa")] != "")==1 & rowSums(tmp[,c("TCRa_1st_cdr_aa","TCRa_2nd_cdr_aa")] != "")==1
#p3 <- ggplot(tsne.plots[specificity.correct$completetable$tcr.frequency >= 10 & tmp.unique.tcr,], aes(x = Coordinate_1, y = Coordinate_2, colour = TCR)) + geom_point(data = tsne.plots[,1:2], colour = "grey", alpha = .2, size = 0.1) + labs(title = "TCRs on tSNE|>=10 cells|donors together") + geom_point(size = 0.3) + facet_wrap(~ TCR) + guides(colour = FALSE) + theme_bw() + theme(strip.text = element_text(size = 5))
#p4 <- ggplot(tsne.plots[specificity.correct$completetable$tcr.frequency >= 10 & tmp.unique.tcr,], aes(x = Coordinate_1, y = Coordinate_2, colour = sampletag)) + scale_color_brewer(palette = "Set1") + geom_point(data = tsne.plots[,1:2], colour = "grey", alpha = .2, size = 0.1) + labs(title = "TCRs on tSNE|>=10 cells|donors seperate") + guides(colour = guide_legend(override.aes = list(size=5))) + geom_point(size = 0.3) + facet_wrap(~ TCR) + theme_bw() + theme(strip.text = element_text(size = 5))
#ggsave("./figures/TCR_tsne_plots.png", marrangeGrob(list(p3,p4),nrow=1,ncol=1),device="png",dpi=300)
#rm(p,p1,p2,p3,p4)

save.image(file = "out/tetramer_analysis.RData")
