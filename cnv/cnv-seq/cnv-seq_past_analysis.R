#!/usr/bin/Rscript
##############################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Thu Jul  4 14:26:54 2013
##############################################
library("cnv")
library("ggbio")

#f<-"jnbc.recal.bam.hits-vs-jnb.recal.bam.hits.log2-0.6.pvalue-0.001.minw-4.count"
data(hg19IdeogramCyto, package = "biovizBase")
cnv.anlysis<-function(bcc.f,chr.len,sample="JHS"){
  bcc.data<-read.delim(bcc.f)
  cnv.summary(bcc.data)
  cnv.print(cnv=bcc.data,file=paste(sample,"_CNV.txt",sep=""))
  #gli1
  #gli1.chr="chr12"
  #gli1.start=57853918
  #gli1.end=57866051
  #plot.cnv(data=bcc.data,chromosome=gli1.chr,from=gli1.start-500000,to=gli1.end+500000)
  
  for(i in c(1:22,"X","Y")){
    chr=paste('chr',i,sep="")
    ideo.p<-plotIdeogram(hg19IdeogramCyto,subchr=chr)
    bcc.p<-plot.cnv(bcc.data,chromosome=paste("chr",i,sep=""),xlabel="")
    
    outfile=paste(sample,"_",chr,".pdf",sep="")
    pdf(file=outfile,width=16,height=10)
    alignPlots(Chr=ideo.p,BCC=bcc.p,heights=c(2,10))
    dev.off()
  }
}
chr.len<-read.table("/Users/wanglab/Google Drive/stanford/wanglab/reference/hg19/hg19.sizes",row.names=1)

cnv.anlysis(bcc.f="RSR.bam.hits-vs-RSN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="RSR")
cnv.anlysis(bcc.f="DLCR1.bam.hits-vs-DLCN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="DLCR1")
cnv.anlysis(bcc.f="DLCR2.bam.hits-vs-DLCN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="DLCR2")
cnv.anlysis(bcc.f="DRR.bam.hits-vs-DRN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="DRR")
cnv.anlysis(bcc.f="RSR.bam.hits-vs-RSN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="RSR")
cnv.anlysis(bcc.f="GRS.bam.hits-vs-GRN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="GRS")
cnv.anlysis(bcc.f="GRR.bam.hits-vs-GRN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="GRR")
#cnv.anlysis(bcc.f="JHS.bam.hits-vs-JHN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="JHS")
cnv.anlysis(bcc.f="JHR.bam.hits-vs-JHN.bam.hits.log2-0.6.pvalue-0.001.minw-4.cnv",chr.len=chr.len,sample="JHR")