#Simple correlation test
#4 Sept 2019
#want to check whether the envrionmental varialbes are correlated

#install pacakages
install.packages("corpcor")
install.packages("ppcor")
install.packages("corrplot")

#load packages
library(corpcor)
library(ppcor)
library(corrplot)

Correlation=read.csv("~/Documents/R_analysis/Fourth corner/FinalMatrices/Site_matrix_27May2019.csv", header=TRUE)
Correlation
Abundance=read.csv("~/Documents/R_analysis/Fourth corner/FinalMatrices/Abundance_matrix_27May2019.csv")
Abundance

#merge specimen and site datasets

BeesCorrelation<-merge(Abundance, Correlation, by=c("X"))

X<-Correlation[,2:10]
cor2pcor(cov(X))

#make correlation as matrix 
Correlation=as.matrix(Correlation)
M=cor(X)

corrplot(M, method = "circle") #plot matrix
write.csv(M,file="CorrelationPlot.csv")
