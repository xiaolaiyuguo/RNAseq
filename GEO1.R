library(org.Hs.eg.db)
g2s<-toTable(org.Hs.egSYMBOL)
head(g2s)
class(g2s)
g2s <- g2s [,-1]
g2s
dim(g2s)

g2e<-toTable(org.Hs.egENSEMBL)
tail(org.Hs.egENSEMBL)
class(g2e)
dim(g2e)

g2 <- cbind(g2s,g2e)
head(g2)

g2$ensembl_id == c(ENSG00000000003.13,
                   ENSG000000000055,
                   ENSG0000000041911,
                   ENSG0000000045712,
                   ENSG0000000046015,
                   ENSG0000000093811)
g2[]

head(g2)


??rbind

install.packages('ALL')









