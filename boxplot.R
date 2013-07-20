########################################
# Boxplot for different size of data
########################################
library(ggplot2)
library(reshape2)

s1<-rnorm(n=10,mean=2,sd=2)  # Normal distribution
s2<-rbinom(n=50,size=10,prob=0.4)  # binomial distribution
s3<-rpois(n=30,lambda=1)     # poisson distribution with mean=1

dat<-data.frame(val=c(s1,s2,s3),dist=c(rep("Normal",length(s1)),
                                       rep("Binominal",length(s2)),
                                       rep("Poisson",length(s3))))

# boxplot
ggplot(data=dat,aes(x=dist,y=val))+geom_boxplot()
