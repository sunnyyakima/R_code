> library(ggfortify)
Loading required package: ggplot2
> 
> mydata_df <- as.data.frame(mydata_log)
> autoplot(prcomp(mydata_df))
> autoplot(prcomp(t(mydata_df)))
> 
