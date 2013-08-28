#!/usr/bin/Rscript
##############################################
# Author: Jiang Li
# Email:  riverlee2008@gmail.com
# Date:   Tue Aug 27 14:10:18 2013
##############################################
loadData<-function(folder){
  files<-list.files(path=folder,pattern="bed$",full.names=TRUE)
  dat<-lapply(files,read.table,sep="\t",stringsAsFactors=FALSE)
  names(dat)<-gsub(pattern=".*\\/(.*?)\\.bed",replacement="\\1",x=files)
  dat
}
dat<-loadData(folder="export_svdzrpkm_svd2/")


cnv.plot<-function(g,extend=1000000){
  # Get data
  sublist<-list()
  for(i in 1:length(dat)){
    d<-dat[[i]]
    start=as.numeric(g[2])-extend
    end=as.numeric(g[3])+extend
    idx<-d$V1==g[1] & d$V2>=start & d$V3<=end
    sublist[[names(dat[i])]]<-d[idx,]
  }
  
  # Plot
  tmp<-sublist[[1]]
  x=(tmp$V2+tmp$V3)/2
  y=tmp$V5
  tt<-smooth.spline(x,y,spar=0.7)
  region=paste(g[1],paste(start,"-",end,sep=""),sep=":")
  plot(x=x,y=y,ylim=c(-3,3),xlim=c(tmp$V2[1],tmp$V3[nrow(tmp)]),xlab="Position",ylab="SVD-ZRPKM values",type="n",pch=1,lty=1,main=region)
  abline(h=0,col='grey',lwd=2)
  for(i in 1:length(sublist)){
    tmp<-sublist[[i]]
    x=(tmp$V2+tmp$V3)/2
    y=tmp$V5
    tt<-smooth.spline(x,y,spar=0.8)
    lines(tt,col=i,lty=i,pch=i)
  }
  legend('topright',legend=names(sublist),col=1:length(sublist),lty=1:length(sublist),bty="n",ncol=5)
  
  # plot gene
  tmp.d<-data.frame()
  count=1
  for(i in 1:nrow(tmp)){
    l<-tmp[i,]
    if(nrow(tmp.d)==0){
      tmp.d<-rbind(tmp.d,l)
    }else{
      if(tmp.d[count,4]==l[,4]){
        tmp.d[count,3]=l[,3]
      }else{
        tmp.d<-rbind(tmp.d,l)
        count=count+1
      }
    }
  }
  tmp.d<-tmp.d[tmp.d$V4!="",]
  
  step=0.2
  
  n1=0
  for(j in 1:nrow(tmp.d)){
    i=j-n1-1 
    if(j%%5==0){
       n1=j/5*5
    }
    rect(xleft=tmp.d[j,2],ybottom=-2-i*step,xright=tmp.d[j,3],ytop=-2-(i-1)*step,col='grey')
    text(x=tmp.d[j,3],y=-2-(i)*(step)+step/2,labels=tmp.d[j,4],cex=0.7,adj=0)
  }
}

gli2<-c("chr2","121493441","121750229")
gli1<-c("chr12","57853918","57866051")
ccnd1<-c("chr11","69455873","69469242")
#ptch1<-c("chr9","98205264","98279247")
extend=100000
pdf("conifier_CNV_gli1_gli2_ptch1_ccnd1.pdf",width=10,height=8)
cnv.plot(gli1,extend=100000)
cnv.plot(gli2)
cnv.plot(ccnd1)
#cnv.plot(ptch1)
dev.off()