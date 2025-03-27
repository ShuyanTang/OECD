#!/usr/bin/Rscript
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(dplyr))

args<-commandArgs(T)

pheno=args[1]


# pheno='ALL'
# pheno='EA'
# pheno='MIX'
# pheno='MA'
# pheno='FF'


input='OE.case_control.LoF_A.matrix'

## ---------- subgroup information ---------------
info = read.delim('allSamples.id',header = T,na.strings='NaN')
if(pheno == 'ALL'){
  cases = info %>% filter(group == 'Case')
}else{
  cases = info %>% filter(group == 'Case') %>% filter(phenotype == pheno)
}

controls=info %>% filter(group == 'Control')

cases.id=gsub("-","\\.",cases$id)
control.id=gsub("-","\\.",controls$id)

case_size=nrow(cases)
ctrl_size=nrow(controls)
size=case_size+ctrl_size


matrix=read.delim(input,header = T, row.names = 1)
data = matrix[,c(control.id,cases.id)]

dat.sum=data.frame(
  ID=rownames(data),
  Control_num=apply(data[,1:ctrl_size],1,sum),
  Case_num=apply(data[,ctrl_size+1:size],1,sum)
)

dat.sum$ID=as.factor(dat.sum$ID)


gid=read.delim('geneids.txt'),header=F)
names(gid)=c('Gene','ID')
gid$ID=as.factor(gid$ID)

d=merge(gid,dat.sum)
d = d %>% filter(Case_num > 0)

## ---------- association ---------------------

for (i in 1:nrow(d)) {
    event.e=d[i,'Case_num']
    n.e=casesize
    event.c=d[i,'Control_num']
    n.c=ctrlsize

    mega=matrix(c(event.e,event.c,n.e-event.e,n.c-event.c),nrow = 2)
    mega.fisher=fisher.test(mega,alternative = 'greater')
    d[i,'fisher.p']=signif(mega.fisher$p.value,3)
    d[i,'fisher.or']=signif(mega.fisher$estimate,3)

}


d=d[order(d$fisher.p),]
write.table(d,paste0('association/asso.',pheno,'.txt'),sep="\t",row.names=FALSE,quote=FALSE)
write.xlsx(d,paste0('association/asso.',pheno,'.xlsx'))




