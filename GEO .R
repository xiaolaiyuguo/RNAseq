BiocManager::install("GEOquery")
library(GEOquery)
browseVignettes("GEOquery")
??GEOquery
help(GEOquery)
gset <- getGEO('GSE42872')
gset.1 <- getGEO('GSE42872',destdir = '.',AnnotGPL =  F ,getGPL = F)
a <- read.table('GSE42872_series_matrix.txt.gz',
                sep = '\t',
                quote = '',
                fill = T,
                comment.char ="!",
                header = T)
str(a)

str(gset.1)

gset.1[[1]]

b = gset.1[[1]]

b

a1<-exprs(b)

class(a1)

rownames(a) <- a[,1]

a <- a[,-1]

str(a)

View(a)

View(a1)

library(hgu95av2.db)

library(org.Hs.eg.db)

b
#GPL6244对应的数据库

BiocManager::install('hugene10sttranscriptcluster.db')
#探针ID和基因ID的转换，基因ID一般用entrez ID
library(hugene10sttranscriptcluster.db)

ids=toTable(hugene10sttranscriptclusterSYMBOL)
head(ids)
table(ids$symbol)
#对对象进行数据统计
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))

exprSet <- a1
table(rownames(exprSet)%in%ids$probe_id)
dim(exprSet)
exprSet<-exprSet[rownames(exprSet)%in%ids$probe_id,]
dim(exprSet)

ids=ids[match(rownames(exprSET),ids$probe_id),]
#match的使用方法
head(ids)
head(exprSet)
exprSet[1:5,1:5]
tmp <- by(exprSet,
          ids$symbol,
          function(x)rownames(x)
          [which.max(rowMeans(x))])
#by函数的使用方法
probes  <- as.character(tmp)
dim(exprSet)
exprSet <- exprSet[rownames(exprSet)%in%probes ,]
dim(exprSet)
exprSet
ids=ids[match(rownames(exprSET),ids$probe_id),]

dim(ids)
dim(exprSet)

#hclust 聚类分析图，用来检测数据
#pca图

tmp<-pData(b)
group_list <- c(rep('control',3),rep('case',3))
group_list


BiocManager::install
library(BiocGenerics)
library(limma)

design <- model.matrix(~0+factor(group_list))
design
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design





