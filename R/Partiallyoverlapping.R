#' The partially overlapping samples t-test
#'
#' Performs a comparison of means using the partially overlapping t-test, for two samples each with paired and unpaired observations.
#' This functions calculates the test statistic, the degrees of freedom, and the p-value. Additionally calculates a confidence interval for the difference in means when requested.
#'  
#' By default, four vectors are to be specified: unpaired observations in Sample 1, unpaired observations in Sample 2, paired observations in Sample 1, paired observations in Sample 2.
#' If the structure of your data is of two vectors, one for each sample, then the option stacked = TRUE can be specified.
#' 
#' @param x1 a vector of unpaired observations in Sample 1 (or all observations in Sample 1 if stacked = "TRUE")
#' @param x2 a vector of unpaired observations in Sample 2 (or all observations in Sample 2 if stacked = "TRUE")
#' @param x3 a vector of paired observations in Sample 1 (not applicable if stacked = "TRUE")
#' @param x4 a vector of paired observations in Sample 2 (not applicable if stacked = "TRUE")
#' @param var.equal a logical variable indicating whether to treat the two variances as being equal. If "TRUE" then the pooled variance is used to estimate the variance, otherwise the Welch approximation to the degrees of freedom is used. 
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param conf.level confidence level of the interval.
#' @param stacked indicator of whether paired and unpaired observations are stacked within one vector ("TRUE"), or if specified as four separate vectors (default)
#' 
#' @return A list which contains the following components of the test:
#' @return statistic  ~~ The value of the t-statistic
#' @return parameter  ~~ The degrees of freedom for the test statistic
#' @return p.value  ~~ The p-value for the test
#' @return estimate ~~ The estimated difference in the means
#' @return conf.int ~~ A confidence interval for the mean appropriate to the specified alternative hypothesis
#' 
#' @details The formula is only applicable for the 2 sample partially overlapping samples t-test. The number of unpaired observations may be zero for up to one of the two samples. The number of paired observations must be of equal length of two or greater. Error messages are generated when these conditions are not true
#' 
#' If your observations are from a datframe of two paired samples with missing observations, use the option "stacked=TRUE". Correpsonding pairs should be given on the same row when this option is applied. 
#' 
#' @examples
#' #The sample means for two groups, "a" and "b" are compared
#' #for a two sided test assuming equal variances.
#' 
#' #Approach 1:
#' #For each sample, unpaired observations and paired observations defined as separate vectors:
#' a.unpaired<-c(20,21,16,18,14,12,14,17)
#' a.paired<-c(14,15,18,20,11,19,14,15)
#' b.unpaired<-c(10,16,18,16,15,14,13,10)
#' b.paired<-c(15,10,15,17,13,19,12,13)
#' Partover.test(a.unpaired,b.unpaired,a.paired,b.paired,var.equal=TRUE) 
#' #p.value = 0.026, the samples from group "a" and group "b" have significantly different means
#'
#' #Equivalently, Approach 2:
#' #Independent observations and the paired samples stacked for each sample:
#' a<-c(20,21,16,18,14,12,14,17,NA,NA,NA,NA,NA,NA,NA,NA,14,15,18,20,11,19,14,15)
#' b<-c(NA,NA,NA,NA,NA,NA,NA,NA,10,16,18,16,15,14,13,10,15,10,15,17,13,19,12,13)
#' Partover.test(a,b,var.equal=TRUE,stacked=TRUE) 
#' #p.value = 0.026, the samples from group "a" and group "b" have significantly different means
#' @export 


Partover.test<-function(x1=NULL,x2=NULL,x3=NULL,x4=NULL,var.equal=FALSE, alternative="two.sided", conf.level=NULL, stacked=FALSE){
  #Seperates observations into format required for stacked = FALSE.
  if (stacked==TRUE){
    if (length(x1)!=length(x2)) stop ("samples must be specified of equal lengths within the matrix for stacked=FALSE. Check structure of data") #The length of a and b must be equal. 
    if ((!is.null(x3))|(!is.null(x4))) stop ("stacked = TRUE option requires only 2 specified vectors, one for each sample. Check structure of data")
    a2<-rep(x1,2)
    b2<-rep(x2,2)
    grp<-rep(NA,(length(x1)*2))
    stacked<-data.frame(a2,b2,grp) 
    for (i in 1:(length (stacked$a2)/2)){if (!is.na(stacked$a2[i])) if (is.na(stacked$b2[i])) stacked$grp[i]<-"aunp"}
    for (i in 1:(length (stacked$a2)/2)){if (!is.na(stacked$a2[i])) if (!is.na(stacked$b2[i])) stacked$grp[i]<-"apaired"}
    for (i in 1:(length (stacked$a2)/2)){if (is.na(stacked$a2[i])) stacked$grp[i]<-"exclude"}  #  The "exclude" coding allows for observations with na on both samples to pass through the system without issue.
    for (i in ((length(stacked$a2)/2)+1): length(stacked$a2)) if (is.na(stacked$a2[i])) if (!is.na(stacked$b2[i])) stacked$grp[i]<-"bunp"
    for (i in ((length(stacked$a2)/2)+1): length(stacked$a2)) if (!is.na(stacked$a[i])) if (!is.na(stacked$b[i])) stacked$grp[i]<-"bpaired"
    for (i in ((length(stacked$a2)/2)+1): length(stacked$a2)) if (is.na(stacked$b[i])) stacked$grp[i]<-"exclude"
    x1<-stacked[stacked$grp=="aunp",1]
    x2<-stacked[stacked$grp=="bunp",2]
    x3<-stacked[stacked$grp=="apaired",1]
    x4<-stacked[stacked$grp=="bpaired",2]}
  
  #Checks that vectors specified are of required length for partially overlapping samples t-test# 
  if (is.null(x1)&&is.null(x2)) {stop ("not enough vectors specified")}
  if (length (x3)!=length(x4)) {stop ("paired observations not of same length")}
  if (length (x3)<2) {stop ("not enough paired observations")}
  
  #Elements of test common to equal and unequal variances assumed
  xbar1 <-mean(c(x1,x3)); xbar2 <-mean(c(x2,x4)); estimate <-xbar1 - xbar2
  n1 <-(length (x1) + length (x3)); n2 <-(length (x2)+length (x4)); n12 <-length(x3)
  if ((stats::sd(x3)==0)|(stats::sd(x4)==0)) {r<-0} else {r<-stats::cor (x3,x4)}
  
  #for equal variances assumed
  if(var.equal==TRUE){
    spooled.1 <- ((n1)-1) *(stats::var(c(x1,x3)))
    spooled.2 <- ((n2)-1) *(stats::var(c(x2,x4)))
    spooled <- sqrt ((spooled.1 + spooled.2) / (n1 +n2 -2))
    denom.1<- (1/n1) +(1/n2)
    denom.2 <- 2*r*n12 / (n1 *n2)
    denom<- spooled *sqrt(denom.1 - denom.2)
    statistic <- estimate / denom
    parameter<-(length (x3) -1) +(((length (x1) + length (x2) + length (x3) -1)/(length (x1) + length (x2) +(2*length (x3))))*(length (x1) + length (x2)))
  }
  
  #for equal variances not assumed  (default)
  if(var.equal==FALSE){
    denom.1<- (stats::var(c(x1,x3))/n1) +(stats::var(c(x2,x4))/n2)
    denom.2 <- (2*r*n12*stats::sd(c(x1,x3))*stats::sd(c(x2,x4))) / (n1 *n2)
    denom<- sqrt(denom.1 - denom.2)
    statistic <- estimate / denom
    welapprox<- (((stats::var(c(x1,x3))/n1)+(stats::var(c(x2,x4))/n2))^2)/((((stats::var(c(x1,x3))/n1)^2)/(n1-1))+(((stats::var(c(x2,x4))/n2)^2 )/(n2-1)))
    parameter <- (n12 -1) + (((welapprox - n12 +1)/(length (x1) + length (x2) +(2*n12)))*( length (x1) + length (x2)))    
  }
  
  #calculates p-value, with p-value fixed in extremes of no variability so that output does not give NaN or NA.
  if (is.na(statistic)) {p.value<-1}
  if (!is.na(statistic)) {if (xor((statistic == -Inf), (statistic == Inf))) p.value<-0 else{
    if(alternative=="less"){
      p.value<-stats::pt(statistic,df=parameter)
    }
    if(alternative=="greater"){
      p.value<-stats::pt(-abs(statistic),df=parameter)
    }
    if(alternative=="two.sided"){
      p.value<-2*stats::pt(-abs(statistic),df=parameter)
    }
  }}
  
  #gives relevant outputs dependent on whether confidence intervals are required
  if(is.null(conf.level)){
    theoutputs<-list(statistic=statistic,parameter=parameter,p.value=p.value,estimate=estimate)
  }
  if(!is.null(conf.level)){
    if(alternative=="two.sided"){
      con<-(1-((1-conf.level)/2))
      critical_val<-stats::qt(con,parameter)
      lower_int<-estimate-(critical_val*denom)
      upper_int<-estimate+(critical_val*denom)
    }
    if(alternative=="less"){
      con<-conf.level +((1-conf.level)/2)
      critical_val<-stats::qt(conf.level,parameter)
      lower_int<-(-Inf)
      upper_int<-estimate+(critical_val*denom)
    }
    if(alternative=="greater"){
      con<-conf.level +((1-conf.level)/2)
      critical_val<-stats::qt(conf.level,parameter)
      lower_int<-estimate-(critical_val*denom)
      upper_int<-Inf
    }
    theoutputs<-list(statistic=statistic,parameter=parameter,p.value=p.value,estimate=estimate,conf.int=c(lower_int,upper_int))
  }
  return(theoutputs) 
}

