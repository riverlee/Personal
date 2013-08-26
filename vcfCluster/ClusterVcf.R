#!/usr/bin/Rscript
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: 2013-08-23 
#############################################
args<-commandArgs(trailingOnly=TRUE)
usage<-function(){
    cat("Usage: ClusterVcf.R genotype.txt cluster.pdf\n")
    q(save="no")
}

if(length(args)!=2){
    usage()
}

infile=args[1]
outfile=args[2]

# Read file
dat<-read.table(infile,sep="\t",comment.char="",header=TRUE)
rownames(dat)<-paste(dat[,1],dat[,2],sep="_")
#only genotype
dat<-dat[,6:ncol(dat)]

#Priciple componet analysis
pdf(file=outfile,width=8,height=8)
pc<-princomp(dat)
load.score<-loadings(pc)
pchs=rep(1,ncol(dat))
cols=pchs
plot(load.score[,1:2],type='n',main="PCA loading plot by princomp")
text(load.score[,1:2],labels=rownames(load.score),pch=pchs,col=cols) # add variable names



#information from http://www.statmethods.net/advstats/factor.html
#By Exploratory Factor Analysis
#The factanal( ) function produces maximum likelihood factor analysis.
mydata<-dat
fit <- factanal(mydata, 3,rotation="varimax")
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n",main="Factor Analysis loading plot by factanal") # set up plot 
text(load,labels=rownames(load),pch=pchs,col=cols) # add variable names
dev.off()

q(save="no")
