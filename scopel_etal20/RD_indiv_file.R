#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).depth", call.=FALSE)
}
depthtab <- fread(args[1], col.names = c("Chromosome", "Position", "Depth"))
c = 0
k = 0
cstart <- 0
cend <- 0
g = sub(".*depth_","",args[1])
pdf(paste0(g,"_SW_SP.pdf"))
med = 2*median(depthtab$Depth, na.rm=TRUE)
a <- unique(depthtab$Chromosome)
a[5:8] <- unique(depthtab$Chromosome)[7:10]
a[9] <- unique(depthtab$Chromosome)[5]
a[10:16] <- unique(depthtab$Chromosome)[11:17]
a[17] <- unique(depthtab$Chromosome)[6]
op <- par(mar=c(1,1,4,1)+0.1, oma = c(5,4,4,1)+0.1, mfcol=c(1,16))
for(i in a){
	if(grepl("chrM", i)){next}
	chrm <- depthtab[depthtab$Chromosome == i,2:3]
  attach(chrm)
  W <- 1:as.integer(max(Position)/1000)*1000
  W <- append(W,max(W)+1000)
  meanDepth <- 0
  for(j in 1:length(W)){meanDepth[j] <- mean(Depth[Position>(W[j]-1000)&Position<=W[j]])}
	par(mar=c(0,0,2,0))
	if(c%%4==0){color="black"}else if(c%%4==1){color="#E44848"}else if(c%%4==2){color="#9FCDF0"}else if(c%%4==3){color="#E0F09F"}
  plot(W[1:max(W)], meanDepth[1:max(W)],pch=19,cex=0.2, ylim=c(0,med), axes = FALSE, col= color)
	title(main = paste0(sub(".*chr", "", i), "\n", median(Depth, na.rm = TRUE), "x"), cex.main = 1)
	k = k + max(W)
	axis(side = 1, at = max(W), labels = k/1000, cex.axis = 1, pos=0)
	if(c==0){axis(side = 2, at = seq(0,med,(med/5)),pos=0)}
  #abline(h=mean(Depth),col="red")
  abline(h=median(Depth, na.rm = TRUE), col="blue")
	abline(h=0, col = "black")
	detach(chrm)
	c=c+1
}
title(main = g, xlab="genomic position (kbp)", ylab="mean depth", outer = TRUE, line = 2)
par(op)
dev.off()
c=0
k=0
rm(list=ls())

