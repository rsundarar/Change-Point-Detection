# Nonparametric Change Point Detection in Multivariate Piecewise
# Stationary Time Series

# Journal of Nonparametric Statistics
# Sundararajan and Pourahmadi (2018+)

# Contact Email: rsundararajan@smu.edu

########################################################################

# List of Required Libraries #

library(boot)
library(np)
library(vars)
library(astsa)
library(MTS)
library(fBasics)
library(cmvnorm)
library(MASS)
library(pracma)
library(gdata)
library(Matrix)

####################################################################################

# Notations and Usage: 

# N = 2^cdl Single value of dyadic neighborhood length choice 

# data1 is the n (x) p input multivariate time series
# p>1 denotes the dimension of the series
# n is the series length

# sig.level: choice of significance leel to be used for the test

####################################################################################

# Use the below two wrapper functions to detect change point locations:

# 1. cpt.locations.multi() - Change pt detection with multiscale choice of N (Section 2.3.1 of paper)
# 2. cpt.locations.single() - Change pt detection with a single choice of N

# Check for the above two function definitions below

#######################################################################################

# Function List Begins #

# L2 Norm
cxnorm=function(x){
  sum1=0
  p=length(x)
  for(i in 1:p){
    sum1=sum1+Mod(x[i])^2  
  }
  return((sum1))
}


##  L2 Metric/Test Statistic - 
# Measure of difference in adjacent spectral matrices
vecdiff=function(data1,time.pt,cdl){
  tt=local.spec(data1,time.pt,cdl)
  
  N = 2^cdl
  p=max(1,dim(data1)[2])
  
  sum1=0
  for(i in 1:N){
    sum1=sum1+
      cxnorm(matrix(tt[[1]][,,i]-tt[[2]][,,i],p*p,1))  
  }
  
  return(sum1/N)
} 

## Finding Local Spectral Matrices - left and right of a time point
local.spec=function(data1,time.pt,cdl){
  
  N=2^cdl
  p=dim(data1)[2]
  m=round(N^(0.5)) # Smoothing
  
  
  data2=data1[(time.pt):(time.pt-N+1),] #Left side of time.pt
  data3=data1[(time.pt+1):(time.pt+N),] #Right side of time.pt
  
  spec.left<-array(0,c(p,p,N))
  spec.right<-array(0,c(p,p,N))
  
  # Uses Daniell kernel #
  ss1=mvspec(data2,kernel("daniell",m),plot=FALSE)$fxx #Left
  ss2=mvspec(data3,kernel("daniell",m),plot=FALSE)$fxx #Right
  
  for(i in 1:N){
    if(i<=N/2) {
      spec.left[,,i] = ss1[,,i]/(2*pi) 
      spec.right[,,i] = ss2[,,i]/(2*pi)
    }
    
    
    if(i>N/2) {
      spec.left[,,i] = t(ss1[,,N-i+1])/(2*pi) 
      spec.right[,,i] = t(ss2[,,N-i+1])/(2*pi)
    }    
    
  }
  
  return(list(spec.left,spec.right))
} # End of function local.spec


# Function to test if a given point is a significant change pt.
# Returns bootstrapped p-value

test.cpt = function(data1,time.pt,cdl){
  
  N = 2^cdl # length of local window
  p = dim(data1)[2]  
  
  tstat.actual = vecdiff(data1,time.pt,cdl) # Actual test stat value
  
  data2=data1[(time.pt-N+1):(time.pt+N),] # Neighborhood of time.pt  
  
  
  #Define Function for test statistic to pass onto tsboot
  vecdiff.stat = function(data22){
    data1[(time.pt-N+1):(time.pt+N),] = data22
    return(vecdiff(data1,time.pt,cdl))  
  }
  
  # Optimal Block length Selection: Politis and White (2004)
  s<-b.star(data2,round=TRUE)[,1]
  # taking the average of the component-wise chosen block length
  l1<-mean(s)
  #Block Bootstrap: Politis and Romano (1992)
  tstat.vector = tsboot(data2,vecdiff.stat,R=500,sim="geom",l=l1)$t
  
  pval = sum(tstat.vector > tstat.actual)/500  
  
  return(pval)    
} # end function test.cpt


#################################################################################################
##################################################################################################
###################################################################################################


# Wrapper functions #

# Function to return change point locations  

# for a single neighborhood length choice N (Needs to be specified)

# Specify 'cdl' that determines N, the dyadic neighborhood length N

cpt.locations.single = function(data1,cdl,sig.level=0.01){
  
  N = 2^cdl # length of dyadic neighborhood length N
  
  p=dim(data1)[2] # Dimension of series 
  
  n1 = dim(data1)[1] #Lenght of series
  
  metric.pts = matrix(0,nrow = n1 , ncol=1)
  cpt.locs = NULL # Final Change Point Locations
  
  for(i in (N+1):(n1-N-1)){
    metric.pts[i,] = vecdiff(data1,i,cdl)
  }
  
  
  flag = 0 # Stopping Indication
  
  while(flag==0){
    
    current.loc = which.max(metric.pts)  # current max location
    
    if(current.loc < N || current.loc > n1-N ) break
    if(metric.pts[current.loc,1]==0) break
    
    current.pval = test.cpt(data1,current.loc,cdl)
    
    if(current.pval<sig.level) {
      cpt.locs = c(cpt.locs,current.loc)
      metric.pts[(current.loc-N+1):(current.loc+N),1] = 0
      next
    }
    
    if(current.pval>sig.level){
      flag=1
      break
    }
    
  } # end of while loop
  
  return(cpt.locs)  
} # end function cpt.locations.single


#  Function to detect change point locations for a multi-scale procedure

# data1 is the n (x) p input multivariate time series

# Adopts a multi-scale procedure for neighborhood length choice N

cpt.locations.multi = function(data1,sig.level=0.01){
  
  p=dim(data1)[2] # Dimension of series  
  n1 = dim(data1)[1] #Lenght of series
  
  N.choices = 2^(1:20)
  
  # Determine set of neighborhood choices N_i,i=1,2,...,n
  lo = sqrt(n1)
  up = n1^(5/6)
  
  # Choices for neighborhood lengths: See Section 2.3.1 of paper
  N.choices = N.choices[which(N.choices >lo & N.choices<up )]
  
  final.cpts = NULL # final list of detected change points
  
  # Initialize array to store change points
  cpts = matrix(0,length(N.choices),1000)
  
  for(k in 1:length(N.choices)){
    cpts.tmp = cpt.locations.single(data1,cdl=log(N.choices[k],base=2))
    if(length(cpts.tmp>0)) cpts[k,1:length(cpts.tmp)] =  cpts.tmp   
  }
  
  
  if(sum(cpts[1,])>0) final.cpts = cpts[1,][which(cpts[1,]>0)] # Initialize
  
  for(k in 2:length(N.choices)){
    
    if(sum(cpts[k,])>0){
      
      l1 = length(which(cpts[k,]>0)) # no.of change pts detected with length N
      
      for(i in 1:l1){
        flag = 0
        
        # Check whether cpts[i] is in the neighborhood of existing change.pts
        if(length(final.cpts)>0){
          for(j in 1:length(final.cpts)){
            nbhd = c(final.cpts[j] - N.choices[k] , final.cpts[j] + N.choices[k])
            if(cpts[k,i] > nbhd[1] & cpts[k,i] < nbhd[2] ) flag = 1
          } # end for over j 
          
          if(flag==0) final.cpts = c(final.cpts , cpts[k,i])
        }
        
        if(length(final.cpts)<0) final.cpts = c(final.cpts,cpts[k,1:l1])  
        
      }  # end for over i
      
    } # end if
    
  } # end for over k
  
  
  return(final.cpts)  
  
} # end function cpt.locations.multi