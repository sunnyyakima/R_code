# PCA plot
library(ggfortify)
mydata_df <- as.data.frame(mydata_log)
autoplot(prcomp(t(mydata_df)))
 
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
