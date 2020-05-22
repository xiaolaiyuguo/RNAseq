rm(list = ls())
#清空当前环境中所有变量
#问题一：
#安装一些R包：
#数据包： ALL, CLL, pasilla, airway 
#软件包：limma，DESeq2，clusterProfiler 
#工具包：reshape2
#绘图包：ggplot2
#不同领域的R包使用频率不一样，在生物信息学领域，尤其需要掌握bioconductor系列包。
install.packages('ALL')
library(ALL)
library(BiocManager)
BiocManager::install('CLL')
BiocManager::install('pasilla')
BiocManager::install('airway')
suppressPackageStartupMessages(library())
BiocManager::install('limma')
BiocManager::install('DESeq2')
BiocManager::install('clusterProfiler')
install.packages('reshape2')
install.packages('ggplot2')

#第二题：了解ExpressionSet对象，比如CLL包里面就有data(sCLLex) ，找到它包含的元素，提取其表达矩阵(使用exprs函数)，查看其大小
#参考：http://www.bio-info-trainee.com/bioconductor_China/software/limma.html
#参考：https://github.com/bioconductor-china/basic/blob/master/ExpressionSet.md
#这个对象其实是对表达矩阵加上样本分组信息的一个封装，由biobase这个包引入。它是eSet这个对象的继承。
#exprs提取表达ExpressionSet的表达矩阵
library(CLL)
data(sCLLex)
sCLLex
class(sCLLex)
str(sCLLex)
exprSet=exprs(sCLLex)
exprSet
class(exprSet)
str(exprSet)
dim(exprSet)
samples=sampleNames(sCLLex)
samples 
#取样本名
pdata=pData(sCLLex)
pdata
#取样本分组
group_list=as.character(pdata[,2])
group_list
#制作样品分组
dim(exprSet)
exprSet[1:5,1:5]



#第三题：了解 str,head,help函数，作用于 第二步提取到的表达矩阵
str(exprSet)
head(exprSet)
help(exprSet)

#第四题：安装并了解 hgu95av2.db 包,看看 ls("package:hgu95av2.db") 后 显示的那些变量
#package:hgu95av2.db包含探针ID的个数据库中基因ID对应的表格
BiocManager::install('hgu95av2.db')
library(hgu95av2.db)
ls("package:hgu95av2.db") 

#第五题：理解 head(toTable(hgu95av2SYMBOL)) 的用法，找到 TP53 基因对应的探针ID
ids=toTable(hgu95av2SYMBOL)
dim(ids)
head(ids)
save(ids,exprSet,pdata,file = 'input.Rdata')
#保存ids，exprSet，pdata，为input.Rdata文件
length(unique(ids$symbol))
head(ids)
tail(ids[,2] == "TP53")
ids[which(ids$symbol=="TP53"),]
ids[ids[,2] == "TP53",]

#第6题：理解探针与基因的对应关系，总共多少个基因，基因最多对应多少个探针，是哪些基因，是不是因为这些基因很长，所以在其上面设计多个探针呢？
length(unique(ids$symbol))
#查看有多少个不重复的symbol
tail(sort(table(ids$symbol)))
#查看每个symbol对应几个探针
table(sort(table(ids$symbol)))
#查看探针重复的次数
plot(table(sort(table(ids$symbol))))

#第7题：第二步提取到的表达矩阵是12625个探针在22个样本的表达量矩阵，找到那些不在 hgu95av2.db 包收录的对应着SYMBOL的探针。
#提示：有1165个探针是没有对应基因名字的。

dim(exprSet)
dim(ids)
head(exprSet)
head(ids)
table(rownames(exprSet) %in% ids$probe_id)
table(ids$probe_id %in% rownames(exprSet))
dim(exprSet)
#exprSet=exprSet[!rownames(exprSet) %in% ids$probe_id,]

#第8题：过滤表达矩阵，删除那1165个没有对应基因名字的探针。
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)

#第9题：整合表达矩阵，多个探针对应一个基因的情况下，只保留在所有样本里面平均表达量最大的那个探针。
#提示，理解 tapply,by,aggregate,split 函数 , 首先对每个基因找到最大表达量的探针。
#然后根据得到探针去过滤原始表达矩阵

ids=ids[match(rownames(exprSET),ids$probe_id),]
head(ids)
exprSet(1:6,1:5)
#match的使用方法
tmp <- by(exprSet,
          ids$symbol,
          function(x)rownames(x)
          [which.max(rowMeans(x))])
#by函数的使用方法
probes  <- as.character(tmp)
dim(exprSet)
exprSet <- exprSet[rownames(exprSet)%in%probes ,]
dim(exprSet)
head(exprSet)
table(table(rownames(exprSet)))

#第10题：把过滤后的表达矩阵更改行名为基因的symbol，因为这个时候探针和基因是一对一关系了。
rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
head(exprSet)

#第11题：对第10步得到的表达矩阵进行探索，先画第一个样本的所有基因的表达量的boxplot,hist,density ， 然后画所有样本的 这些图
#A:参考：http://bio-info-trainee.com/tmp/basic_visualization_for_expression_matrix.html
#B:理解ggplot2的绘图语法，数据和图形元素的映射关系
#boxplot，hist，denstity
exprSet[1:5,1:5]
library(reshape2)
exprSet_L <- melt(exprSet)
head(exprSet_L)
table(exprSet_L$Var1 == 'MAPK7')
colnames(exprSet_L)=c('symbol','sample','value')
head(exprSet_L)
exprSet_Z<-exprSet_L
exprSet_Z$group=rep(group_list,each=nrow(exprSet))
#rep函数的使用
head(exprSet_Z,20)

#boxplot
library(ggplot2)
p=ggplot(exprSet_Z[exprSet_Z$sample=='CLL11.CEL',],aes(x=sample,y=value,fill=group))+geom_boxplot()
print(p)
#hist
p=ggplot(exprSet_Z[exprSet_Z$sample=='CLL11.CEL',],aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
print(p)
#denstity
p=ggplot(exprSet_Z[exprSet_Z$sample=='CLL11.CEL',],aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
print(p)

#第12题：理解统计学指标mean,median,max,min,sd,var,mad并计算出每个基因在所有样本的这些统计学指标，最后按照mad值排序，取top 50 mad值的基因，得到列表。
#注意：这个题目出的并不合规，请仔细看。
#R语言中的apply函数
??apply
#函数定义：apply(X, MARGIN, FUN, ...)
#参数列表：X:数组、矩阵、数据框
#MARGIN: 按行计算或按按列计算，1表示按行，2表示按列
#FUN: 自定义的调用函数
#…: 更多参数，可选
g_mean <- tail(sort(apply(exprSet,1,mean)),50)
g_median <- tail(sort(apply(exprSet,1,median)),50)
g_max <- tail(sort(apply(exprSet,1,max)),50)
g_min <- tail(sort(apply(exprSet,1,min)),50)
g_sd <- tail(sort(apply(exprSet,1,sd)),50)
g_var <- tail(sort(apply(exprSet,1,var)),50)
g_mad <- tail(sort(apply(exprSet,1,mad)),50)
g_mad
names(g_mad)
#mad所有数减去中位数的绝对值，取中位数，和标准差相对应

#13题：根据第12步骤得到top 50 mad值的基因列表来取表达矩阵的子集，并且热图可视化子表达矩阵。试试看其它5种热图的包的不同效果。
help(pheatmap)
??pheatmap
library(pheatmap)
choose_gene=names(tail(sort(apply(exprSet,1,mad)),50))
choose_matrix=exprSet[choose_gene,]
choose_matrix
choose_matrix[1:5,1:5]
choose_matrix_1 <- t(choose_matrix)
#行列转换
choose_matrix_1[1:5,1:5]
choose_matrix_2=scale(t(choose_matrix))
#R语言中scale函数的意义
#scale函数是将一组数进行处理，默认情况下是将一组数的每个数都减去这组数的平均值后再除以这组数的均方根。
choose_matrix_2[1:5,1:5]
choose_matrix_3=t(scale(t(choose_matrix)))
choose_matrix_3[1:5,1:5]
install.packages('pheatmap')
pheatmap(choose_matrix_3)


#第14题：取不同统计学指标mean,median,max,mean,sd,var,mad的各top50基因列表，使用UpSetR包来看他们之间的overlap情况。
install.packages('UpSetR')
library(UpSetR)
g_all <- unique(c(names(g_mean),names(g_median),names(g_max),names(g_min),
                  names(g_sd),names(g_var),names(g_mad) ))
dat=data.frame(g_all=g_all,
               g_mean=ifelse(g_all %in%  names(g_mean) ,1,0),
               g_median=ifelse(g_all %in%  names(g_median) ,1,0),
               g_max=ifelse(g_all %in%  names(g_max) ,1,0),
               g_min=ifelse(g_all %in%  names(g_min) ,1,0),
               g_sd=ifelse(g_all %in%  names(g_sd) ,1,0),
               g_var=ifelse(g_all %in%  names(g_var) ,1,0),
               g_mad=ifelse(g_all %in%  names(g_mad) ,1,0))
upset(dat,nsets = 7)

#第15题：在第二步的基础上面提取CLL包里面的data(sCLLex) 数据对象的样本的表型数据。
pdata=pData(sCLLex)
group_list=as.character(pdata[,2])
group_list
dim(exprSet)
exprSet[1:5,1:5]

#第16题：对所有样本的表达矩阵进行聚类并且绘图，然后添加样本的临床表型数据信息(更改样本名)
colnames(exprSet)=paste(group_list,1:22,sep='')
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
                cex = 0.7, col = "blue")
hc=hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10)) 
plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)

#第17题：对所有样本的表达矩阵进行PCA分析并且绘图，同样要添加表型信息。
install.packages('ggfortify')
library(ggfortify)
exprSet[1:5,1:5]
df=as.data.frame(t(exprSet))
df[1:5,1:5]
head(df)
df$group=group_list 
df$group
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')

install.packages('FactoMineR')#画主成分分析图需要加载这两个包
install.packages('factoextra') 
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") 
df=as.data.frame(t(exprSet))
dat.pca <- PCA(df, graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group_list, # color by groups
             # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)


#第18题：根据表达矩阵及样本分组信息进行批量T检验，得到检验结果表格
dat = exprSet
group_list=as.factor(group_list)
group1 = which(group_list == levels(group_list)[1])
group2 = which(group_list == levels(group_list)[2])
dat1 = dat[, group1]
dat2 = dat[, group2]
dat = cbind(dat1, dat2)
pvals = apply(exprSet, 1, function(x){
  t.test(as.numeric(x)~group_list)$p.value
})
p.adj = p.adjust(pvals, method = "BH")
avg_1 = rowMeans(dat1)
avg_2 = rowMeans(dat2)
log2FC = avg_2-avg_1
DEG_t.test = cbind(avg_1, avg_2, log2FC, pvals, p.adj)
DEG_t.test=DEG_t.test[order(DEG_t.test[,4]),]
DEG_t.test=as.data.frame(DEG_t.test)
head(DEG_t.test)
#第19题：使用limma包对表达矩阵及样本分组信息进行差异分析，得到差异分析表格，重点看logFC和P值，画个火山图(就是logFC和-log10(P值)的散点图。)。
suppressMessages(library(limma)) 
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design

## 下面的 contrast.matrix 矩阵非常重要，制定了谁比谁这个规则
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix 
##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
## volcano plot
DEG=nrDEG
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
head(DEG)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)

#第20题：对T检验结果的P值和limma包差异分析的P值画散点图，看看哪些基因相差很大。
### different P values 
head(nrDEG)
head(DEG_t.test)
DEG_t.test=DEG_t.test[rownames(nrDEG),]
plot(DEG_t.test[,3],nrDEG[,1]) ## 可以看到logFC是相反的
plot(DEG_t.test[,4],nrDEG[,4]) # 可以看到使用limma包和t.test本身的p值差异尚可接受
plot(-log10(DEG_t.test[,4]),-log10(nrDEG[,4]))

exprSet['GAPDH',]
exprSet['ACTB',]
exprSet['DLEU1',]


library(ggplot2)
library(ggpubr)
my_comparisons <- list(
  c("stable", "progres.")
)
dat=data.frame(group=group_list,
               sampleID= names(exprSet['DLEU1',]),
               values= as.numeric(exprSet['DLEU1',]))
ggboxplot(
  dat, x = "group", y = "values",
  color = "group",
  add = "jitter"
)+
  stat_compare_means(comparisons = my_comparisons, method = "t.test")


## heatmap 
library(pheatmap)
choose_gene=head(rownames(nrDEG),25)
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)



