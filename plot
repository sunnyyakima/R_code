# boxplot with points
data1 <- read.csv("corr_matrix_4K_afterfilter_t.txt", header=T, row.names=1)
data2 <- read.csv("matched_index.csv", header=F, row.names=1) 
pdf("age_boxplot_with_corr.pdf")
boxplot(data1)
points(1:ncol(data1),do.call(c,lapply(data2,function(x){c(data1[x[2]-2,x[1]])})),col="red",pch=16)
dev.off()

# scatter plot
y <- read.csv("Y_genes3_t.csv", row.names=1)
library(ggplot2)
qplot(y$RPS4Y1,y$DDX3Y,col=y$X.gender.)
y <- read.csv("IG_expr_probes2_t.csv", header=T, row.names=1)
pairs(y) # pair-wise plot


# PCA plot
library(ggfortify)
mydata_df <- as.data.frame(mydata_log)
autoplot(prcomp(t(mydata_df)))

data <- read.csv("ExprM_w_clinical_t.csv", header=T, row.names=1, colClasses=c(rep("character",12), rep("numeric",33781)))
expr <- data[,-(1:11)]
expr <- as.data.frame(expr)
autoplot(prcomp(expr))
df <- as.data.frame(data)
autoplot(prcomp(expr), data = df, colour = 'site')


# scatter plot
library(ggplot2)
ggplot(bcl2l_df, aes(x=x, y=bcl2_score) ) + geom_point(size=2, shape=23)
bcl2l_df$Cancer_type <- as.factor(bcl2l_df$Cancer_type)
ggplot(bcl2l_df, aes(x=x, y=bcl2_score, color=Cancer_type) ) + geom_point(size=2, shape=23)
 
# correlation
library("ggpubr")
mydata <- read.csv("cmp_predict_tim.csv", header=T, row.names=1)
ggscatter(mydata, x="tim_result", y="predict", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", xlab = "tim_crispr_results", ylab= "glmnet predict")
pdf("cmp_glmnet_tim.pdf")
ggscatter(mydata, x="tim_result", y="predict", add = "reg.line", cor.coef = TRUE, cor.method = "pearson", xlab = "tim_crispr_results", ylab= "glmnet predict")
dev.off()

# heatmap
setwd("/Users/yuliu/Bladder/publicData/heatmap3")
library("gplots")
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
 
data_in <- read.csv("public_EMT_data3.csv", header=T, row.names=1)
data_matrix <- data.matrix(data_in)
 
x <- scale(t(as.matrix(data_matrix)))
x <- t(x)
 
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
 
grp = c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue")
grp = c("cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta")
subclass=c("yellow","yellow","yellow","red","red","yellow","yellow","red","blue","blue","red","yellow","blue","red","yellow","red","yellow","red","red","blue","yellow","yellow","yellow","red","red","yellow","yellow","red","blue","blue","red","yellow","blue","red","yellow","red","yellow","red","red","blue")
subclass=c("grey","grey","grey","darkblue","darkblue","grey","grey","darkblue","green","green","darkblue","grey","green","darkblue","grey","darkblue","grey","darkblue","darkblue","green","grey","grey","grey","darkblue","darkblue","grey","grey","darkblue","green","green","darkblue","grey","green","darkblue","grey","darkblue","grey","darkblue","darkblue","green") 

clab = cbind(grp, subclass)
colnames(clab)=c("condition","subclass")

heatmap.3(x, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12), Rowv=TRUE, Colv=TRUE, ColSideColors=clab, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", labCol=FALSE, col=rev(heat.colors(75)), ColSideColorsSize=2)

library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
>  
data_in <- read.csv("public_whole_set_log2_FC.csv", header=T, row.names=1)
data_matrix <- data.matrix(data_in)
column_annotation <- c("grey","grey","grey","darkblue","darkblue","grey","grey","darkblue","green","green","darkblue","grey","green","darkblue","grey","darkblue","grey","darkblue","darkblue","green")
column_annotation <- as.matrix(column_annotation)
colnames(column_annotation) <- c("subclass")

pdf("EMT_pathway_log_FC_whole2.pdf")
heatmap.3(data_matrix, ColSideColors=column_annotation)
dev.off()

# heatmap.2
setwd("/Users/yuliu/Bladder/publicData_GSE87304/Sinai_subtype/neuronal_heatmap")
 mydata <- read.csv("neuronal_topgenes.csv", header=T, row.names=1)
library(gplots)
library(RColorBrewer)
data_matrix <- data.matrix(mydata)
x <- scale(t(as.matrix(data_matrix)))
x <- t(x)
heatmap.2(x, scale = "none", trace = "none", col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), distfun = function(x) as.dist(1-cor(t(x))), hclustfun = function(x) hclust(x, method="ave"), margins=c(4,7), Colv=FALSE)

# ordered gene list
h <- heatmap.2(x, scale = "none", trace = "none", col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), distfun = function(x) as.dist(1-cor(t(x))), hclustfun = function(x) hclust(x, method="ave"), margins=c(4,7), Colv=TRUE)
rownames(x)[h$rowInd]


# stacked histogram
exprM <- read.csv("IGH_genes2_t.csv", header=T)
igh_class <- read.csv("GSE26863_IGH_class.csv", header=T)
merged <- merge(exprM, igh_class, by="MMRCid")
df <- as.data.frame(merged)
library(ggplot2)
qplot(IGHM, data=df, fill=HeavyChainClass, binwidth=0.2)

igh_class$HeavyChainClass[igh_class$HeavyChainClass != "IgG"] = "nonIgG"
merged <- merge(exprM, igh_class, by="MMRCid")
df <- as.data.frame(merged)
qplot(IGHG1, data=df, fill=HeavyChainClass, binwidth=0.2)

# more histgram
http://www.sthda.com/english/wiki/ggplot2-histogram-easy-histogram-graph-with-ggplot2-r-package

igh_class <- read.csv("GSE26863_IGH_class.csv", header=T)
igh_class$HeavyChainClass[igh_class$HeavyChainClass != "IgG"] = "nonIgG"  # not ideal
head(igh_class)
    MMRCid HeavyChainClass
1 MMRC0001            <NA>
2 MMRC0002             IgG
3 MMRC0003             IgG
4 MMRC0006            <NA>
5 MMRC0009             IgG
6 MMRC0010             IgG
> 
ggplot2.histogram(data=df, xName='IGHG1', binwidth=0.2, groupName='HeavyChainClass', legendPosition="top",alpha=0.5, addDensity=TRUE)

igh_class <- read.csv("GSE26863_IGH_class.csv", header=T)
head(igh_class)
sele = igh_class$HeavyChainClass != 'IgG'
igh_class[sele,]
sele [is.na(sele)] = T
igh_class[sele,]
tmp = c('isIgG', 'notIgG')[as.numeric(sele)+1]
igh_class$groupIgG = tmp
head(igh_class)

