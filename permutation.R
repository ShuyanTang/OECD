#!/usr/bin/Rscript
args<-commandArgs(T)

pheno=args[1]
n=args[2]

input='OE.case_control.LoF_A.matrix'
output=paste('permutation/perm',n,pheno,sep=".")

## ---------- subgroup information ---------------
info = read.delim('allSamples.id',header = T,na.strings='NaN')
if(pheno == 'all'){
  cases = info %>% filter(group == 'Case')
}else{
  cases = info %>% filter(group == 'Case') %>% filter(phenotype == pheno)
}

controls=info %>% filter(group == 'Control')

cases.id=gsub("-","\\.",cases$id)
control.id=gsub("-","\\.",controls$id)
all.id=c(cases.id,controls.id)

case_size=nrow(cases)
ctrl_size=nrow(controls)
size=case_size+ctrl_size


md<-read.delim(input,header = T)
data = md[,names(md) %in% all.id]


## ---------- permutation ---------------------

dir="permutation"
setwd(dir)

permdf<-data.frame()
for (i in 1:nrow(data)){
    d=data[i,1:(size)]
    PermSamples <- sample(d, size= size, replace=FALSE)
    
    exp.ctrl<-sum(PermSamples[(case_size+1):size])
    exp.case<-sum(PermSamples[1:case_size])
    m<-matrix(c(exp.case,exp.ctrl,case_size-exp.case,ctrl_size-exp.ctrl),
                nrow=2,ncol=2)
    f<-fisher.test(m,alternative = 'greater')
    permdf[i,'ID']<-md[i,'ID']
    permdf[i,'P']<-f$p.value
}

library(dplyr)
permdf=permdf %>% arrange(P)

write.table(permdf,output,row.names=F,col.names=F,sep="\t")

