#!/usr/bin/Rscript
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#############################################
# Take a gtf file and a list of bed 
# files to calculate the count data and rpkm
args<-commandArgs(trailingOnly=TRUE)

############################
# parameter checking
############################
if(length(args)<2){
  cat("Usage: bed2count.R gene.gtf in.bed [in2.bed] [in3.bed]\n")
  quit(save="no",status=1)
}

gtf<-args[1]
args<-args[-1]

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(GenomicRanges))
############################
# Functions
############################
# 1) Fun1, Read gtf file
gtf2GRangesList <- function(myfile=null) {
  gtf <- read.table(myfile, header=FALSE,sep="\t",stringsAsFactors=FALSE)
  colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", 
                     "strand", "frame","attributes")
  # change according to the genome
  chronly <- c(paste("chr",1:22,sep=""), "chrX", "chrY")
  # Cleanup to remove non-chromosome rows
  gtf <- gtf[as.character(gtf$seqname) %in% chronly, ] 
  
  # only get the exon region
  gtf<-gtf[gtf$feature == "exon",]
  # read the attributes columns, here gene_id and gene_name have the save values,
  # get gene_id from attributes column
  gene.ids <-  gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes) 
  # get gene_name from attributes column
  # gene.name <- gsub(".*gene_name (.*?);.*", "\\1", gtf$attributes) 
  # get transcript_id from attributes column
  transcript.ids<- gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes) 
  # skip those have no value
  gene.index<-gene.ids!=""                    
  
  gene.gr<-GRanges(seqnames=gtf$seqname[gene.index],
                   ranges=IRanges(gtf$start[gene.index],gtf$end[gene.index]),
                   #strand=gtf$strand[gene.index],
                   tx_id=transcript.ids[gene.index],
                   gene_id=gene.ids[gene.index])
  gene.gr.list<-split(gene.gr,gene.ids[gene.index])
  # collapse exon regions for genes having multipl transcripts/isoform
  # gene.gr.list<-lapply(gene.gr.list,reduce)
  return(gene.gr.list)
}

# Fun2: get count based on bam oject  and grangelist
getCounts<-function(aln,grangelist){
  counts<-countOverlaps(grangelist,aln)
  names(counts)<-names(grangelist)
  return(counts)
}

#######################
# Main program
#######################
cat("[",date(),"] Loading transcriptom annotation ...\n")
grangelist<-gtf2GRangesList(gtf)
gene.exon.length<-sum(width(reduce(grangelist)))

for(bed in args){
  n<-basename(bed)
  count.out<-paste(n,"_count.txt",sep="")
  rpkm.out<-paste(n,"_rpkm.txt",sep="")
  if(file.exists(bed)){
    cat("[",date(),"] Doing bed ...\n")
    con<-NULL
    if(grepl("gz$",x=bed,perl=TRUE)){
      con<-gzfile(bed)
    }else{
      con<-file(bed)
    }
    
    aln <- import.bed(con)
    exp <- getCounts(aln,grangelist)
    write.table(exp,file=count.out,sep="\t",col.names=FALSE,quote=FALSE)
    #RPKM
    rpkm<-exp/gene.exon.length
    rpkm<-rpkm/sum(exp)
    rpkm<-rpkm*1e9
    write.table(rpkm,file=rpkm.out,sep="\t",col.names=FALSE,quote=FALSE)
    rm(aln,exp,rpkm)
    
  }else{
    cat("[",date(),"] Skip",bed,", it not exists\n")
  }
}
