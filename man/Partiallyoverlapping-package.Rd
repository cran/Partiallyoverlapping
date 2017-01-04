\name{Partiallyoverlapping-package}
\alias{Partiallyoverlapping-package}
\alias{Partiallyoverlapping}
\docType{package}
\title{
Partially overlapping samples t-tests  
~~ Partiallyoverlapping ~~
}
\description{
The "partially overlapping samples t-tests", for the comparison of means for two samples which include both paired observations and independent observations, can be performed.
}
\details{
\tabular{ll}{
Package: \tab Partiallyoverlapping\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-01-01\cr
License: \tab GPL-3 \cr
}
~~ This package should be used when comparing the means of two samples, that include at least 2 pairs, and a number of unpaired observations~~
}
\author{
Ben Derrick <ben.derrick@uwe.ac.uk>

Maintainer: Ben Derrick <ben.derrick@uwe.ac.uk>
}
\references{
Derrick, B., Russ, B., Toher, D. & White P (2017). Test statistics for the comparison of means for two samples which include both paired observations and independent observations. Journal of Modern Applied Statistical Methods, 16(1).
}
\examples{
#The sample means for two groups, "a" and "b" are compared
#for a two sided test assuming equal variances.
#
#Approach 1:
#For each sample, unpaired observations and paired observations defined as separate vectors:
a.unpaired<-c(20,21,16,18,14,12,14,17)
a.paired<-c(14,15,18,20,11,19,14,15)
b.unpaired<-c(10,16,18,16,15,14,13,10)
b.paired<-c(15,10,15,17,13,19,12,13)
Partover.test(a.unpaired,b.unpaired,a.paired,b.paired,var.equal=TRUE) 
#p.value = 0.026, the samples from group "a" and group "b" have significantly different means
#
#Equivalently, Approach 2:
#Independent observations and the paired samples stacked for each sample:
a<-c(20,21,16,18,14,12,14,17,NA,NA,NA,NA,NA,NA,NA,NA,14,15,18,20,11,19,14,15)
b<-c(NA,NA,NA,NA,NA,NA,NA,NA,10,16,18,16,15,14,13,10,15,10,15,17,13,19,12,13)
Partover.test(a,b,var.equal=TRUE,stacked=TRUE) 
#p.value = 0.026, the samples from group "a" and group "b" have significantly different means
}
