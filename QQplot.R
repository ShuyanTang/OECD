#!/usr/bin/Rscript
library(ggplot2)
library(ggrepel)
library(ggsci)
library(dplyr)


args<-commandArgs(T)

pheno=args[1]


## ---------- read association data ---------------
asso<-read.delim(paste0('association/asso.',pheno,'.txt'),header = T)
# names(df)=c('Gene','ID','Case_num','Control_num','P','OR')

## ---------- read permutation data --------------- 
perm1=read.delim(paste0('permutation/perm.1.',pheno),header=F)
perm1 = perm1 %>% arrange(V2)

perm.data=data.frame(perm1$V2)
for (i in 2:1000){
  df=read.delim(paste('permutation/perm',n,pheno,sep="."),header = F)
  tmp=df %>% arrange(V2)
  perm.data=cbind(perm.data,tmp$V2)
}

write.table(perm.data,paste0('permutation/perm1000.summary.',pheno))

qq.data = data.frame(
    gene = asso$Gene,
    id = asso$ID,
    observed = asso$P,
    expected = -log10(sort(apply(perm.data,1,mean))) 
)


## ---------- QQ plot ---------------------------  
y_threshold_value = -log10(0.001)
qq.plot <- ggplot(qq.data) +
  geom_point(aes(expected, observed, shape = factor(shape), size = size, color = color)) +
  scale_shape_manual(values = c("16" = 16, "17" = 17)) +
  scale_size_identity() +
  scale_shape_manual(values = c(16, 17), labels = c("Known", "Novel")) +
  scale_color_manual(values = c("#9F20F0" = "#9F20F0", "#D95F0E" = "#D95F0E", "grey" = "grey")) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, color = "red") +
  geom_hline(yintercept = y_threshold_value, linetype = "dashed", color = '#f4c0a6') +
  geom_hline(yintercept = -log10(0.05/nrow(qq.data)), linetype = "dashed", color = "#D95F0E") +
  geom_text_repel(data = qq.data %>% filter(observed > y_threshold_value), 
                  aes(x = expected, y = observed, label = Gene), 
                  hjust = 0, nudge_x = -2, nudge_y = 0.5, max.overlaps = 20) +
  xlab(expression(-log[10](Expected ~ p))) +
  ylab(expression(-log[10](Observed ~ p))) +
  scale_x_continuous(limits = c(0, 4), expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0, 0.04)),
                     breaks = seq(0, 10, by = 2), 
                     labels = c(seq(0, 8, by = 2), ">10")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y.left = element_text(color = "black"),
        axis.line.y.left = element_line(color = "black"),
        axis.line.x.bottom = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None") +
  ggtitle(pheno)

ggsave(paste0('plot/qqplot.',pheno,'.pdf'),qq.plot,width=6,height=6,useDingbats=FALSE)


