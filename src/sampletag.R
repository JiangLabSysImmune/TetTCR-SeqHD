library(mixdist)
library(ggplot2)

args <- commandArgs(TRUE)
sampletag <- read.csv(args[1], row.names=1)
rhapsody <- read.csv(args[2], row.names=1, comment.char="#")
miss <- setdiff(rownames(rhapsody),rownames(sampletag))
sampletag[miss,] <- 0
sampletag <- sampletag[rownames(rhapsody),]
sampletag.max <- apply(sampletag,1,which.max)
sampletag.det <- unique(sampletag.max)
names(sampletag.det) <- colnames(sampletag)[sampletag.det]
sampletag.log <- log10(sampletag+1)
thr <- numeric(length = length(sampletag.det))
names(thr) <- colnames(sampletag)[sampletag.det]

set.seed(0)
pdf(file="./figures/sampletag_thresholds.pdf",width = 8,height=8)
layout(matrix(c(1:4),2,2,byrow = T))
for(i in 1:length(sampletag.det)){
  his <- hist(sampletag.log[,sampletag.det[i]], breaks=30, plot=FALSE)
  df <- data.frame(mid=his$mids,cou=his$counts)
  for(j in seq(0.1,5,by=0.1)){
    print(j)
    res <- try(mix.model <- mix(as.mixdata(df),mixparam(mu=c(0,2),sigma=j)))
    if(inherits(res, "try-error")){
      next
    }else{
      break
    }
  }
  plot(res,main=paste(names(sampletag.det)[i]))
  f <- function(x){dnorm(x,res$parameters$mu[1],res$parameters$sigma[1]) * res$parameters$pi[1] -dnorm(x,res$parameters$mu[2],res$parameters$sigma[2]) * res$parameters$pi[2]}
  thr[i] <- uniroot(f,interval = c(res$parameters$mu[1],res$parameters$mu[2]))$root
  abline(v=thr[i], col="blue")
}
dev.off()
print(thr)

samples <- list()
for(i in 1:length(sampletag.det)){
  samples[[i]] <- rownames(sampletag)[sampletag.log[,sampletag.det[i]] > thr[i]]
}
names(samples) <- names(thr)
print(summary(samples))

sampletag.call <- data.frame(row.names = rownames(rhapsody))
sampletag.call$call <- sapply(rownames(sampletag.call), function(i) {
  paste(names(unlist(lapply(samples, function(x) grep(paste("^",i,"$",sep=""), x)))),collapse = "|")
  })
emp <- which(sampletag.call$call=="")
sampletag.call$call[emp] <- sapply(emp, function(i) {paste(colnames(sampletag)[which(sampletag[i,] > 0.5*sum(sampletag[i,]))],collapse = "|")})

sampletag.call$tag.final <- sapply(sampletag.call$call, function(i) {
  ifelse(grepl("st1$",i) | grepl("str1\\|",i), return("tet_neg"),
         ifelse(grepl("\\|",i), return("multiplets"), ifelse(i=="", return("undetermined"), return(i))))
  })


sampletag.pca <- prcomp(sampletag.log)
pca.res <- as.data.frame(sampletag.pca$x)
ggplot(pca.res,aes(PC1,PC2))+geom_point(aes(color=as.factor(sampletag.call$tag.final)),size=1)+labs(colour = "SampleTags")+coord_fixed()+ theme_bw()
ggsave("./figures/sampletag_pca.pdf",width=5,height=4)
write.csv(sampletag.call,"./out/sampletag.csv")
