#' Calculates the value of the shrinkage coefficient for the off-diagonal matrices
#'
#' @description Calculates the value of the shrinkage coefficient for the off-diagonal
#' matrices as described in \[Cruz-Cano et al., 2012\]
#'
#' @param xs Matrix with the values for the datasets X and Y.
#' @param p Number of variables in the dataset X.
#' @param q Number of variables in the dataset Y.
#'
#' @return Shrinkage coefficient for the off-diagonal matrices used to calculate the FRCC canonical structure.
#' @references Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis,
#'             Computational Statistics & Data Analysis, Volume 70, 2014, Pages 88-100,
#'             ISSN 0167-9473, https://doi.org/10.1016/j.csda.2013.09.020.
#'
#' @keywords Internal
#' @author Raul Cruz-Cano
#' @export
off.diagonal.lambda <- function(xs, p,q)
  #xs must be standarized data, which in our case is performed in frcc but just in case
{  xs<-as.matrix(xs) # standardize input matrix by standard deviations
# xs = wt.scale(xs, w, center=TRUE, scale=TRUE) # standardize data matrix   #comes from wt.scale.R
#at this point xs is cbind(X,Y)
w=matrix(data=c(1),ncol=dim(xs)[1],nrow=1)#all experimental units will be weigthed equally
# bias correction factors
w2 = sum(w*w)           # for w=1/n this equals 1/n   where n=dim(xs)[1]
h1w2 = w2/(1-w2)        # for w=1/n this equals 1/(n-1)

sw = sqrt(w)
Q1.squared = (crossprod(sweep(xs, MARGIN=1, STATS=sw, FUN="*")))^2
Q1.squared_XY = Q1.squared[1:p,p+1:q]#extrating part of XY
Q2 = crossprod(sweep(xs^2, MARGIN=1, STATS=sw, FUN="*")) - Q1.squared
Q2_XY = Q2[1:p,p+1:q]#extrating part of XY
#removing the part about substracting the diagonal because now the target is the null-matrix, not the identity
denominator = sum(Q1.squared_XY)#-sum(diag(Q1.squared))
numerator = sum(Q2_XY)#-sum(diag(Q2))

if(denominator == 0)
  lambda = 1
else
  lambda = min(1, numerator/denominator * h1w2)

return (lambda)
}
#====================function needed to rearrange the res.frcc structure if th correlations are out of order due to insufficient data ===========================
#' Rearranges the canonical structure according to the canonical correlations
#'
#' By using the minimum risk estimators of the correlation matrices instead of the
#'   sample correlation matrices the FRCC algorithm might disrupt the order of the
#'   canonical correlations and hence of the canonical structure. This is unacceptable
#'   for the algorithm used  to calculate the p-values which requires the canonical
#'   correlations to be ordered in a descending order. This function rearranges the
#'   canonical structure according to the canonical correlations from largest to smallest.
#'
#' @param res.frcc List containing a canonical structure produced by the function frcc.
#'
#' @return A list containing the sorted canonical structure.
#'
#' @references Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis,
#'             Computational Statistics & Data Analysis, Volume 70, 2014, Pages 88-100,
#'             ISSN 0167-9473, https://doi.org/10.1016/j.csda.2013.09.020.
#'
#' @author Raul Cruz-Cano
#' @keywords Internal
#' @export
rearrange.frcc <-function(res.frcc)
{#  print(res.frcc)
  #Dimensions of the problem
  p<-dim(res.frcc$canonical_weights_X)[1]
  #p = #number of variables in X
  q<-dim(res.frcc$canonical_weights_Y)[1]
  #q = #number of variables in Y
  aux.res.frcc<-list(cor=0,canonical_weights_X=matrix(data=c(0),nrow=p,ncol=1),canonical_weights_Y=matrix(data=c(0),nrow=q,ncol=1))
  for (i in 1:(q-1))
  {
    for (j in (i+1):q)
    {  #we have found a correlation grater than correlation i
      if (res.frcc$cor[j] > res.frcc$cor[i])
      { #print("#swaping the correlations and the corresponding coefficients (at this point we don't have p-values or canonical factor loading to worry about")
        aux.res.frcc$cor <-  res.frcc$cor[i]
        aux.res.frcc$canonical_weights_X <-  res.frcc$canonical_weights_X[,i]
        aux.res.frcc$canonical_weights_Y <-  res.frcc$canonical_weights_Y[,i]

        res.frcc$cor[i]  <-  res.frcc$cor[j]
        res.frcc$canonical_weights_X[,i] <- res.frcc$canonical_weights_X[,j]
        res.frcc$canonical_weights_Y[,i] <- res.frcc$canonical_weights_Y[,j]

        res.frcc$cor[j] <- aux.res.frcc$cor
        res.frcc$canonical_weights_X[,j] <- aux.res.frcc$canonical_weights_X
        res.frcc$canonical_weights_Y[,j] <- aux.res.frcc$canonical_weights_Y
      }#end if
    }#end for j now the CCA i is in the correct place
  }#end for i
  res.frcc
} #end rearrange.frcc function

#=========================function which performs the frccA ======================================================================

#' This function implements the Fast Regularized Canonical Correlation Analysis
#'
#' @description This function implements the Fast Regularized Canonical Correlation algorithm
#' described in \[Cruz-Cano et al., 2014\].
#'
#' The main idea of the algorithm is using the minimum risk estimators of the correlation
#'  matrices described in \[Schafer and Strimmer, 2008\] during the calculation of the Canonical
#'  correlation Structure.
#'
#' It can be considered an extension of the work for two set of variables (blocks)
#' mentioned in \[Tenenhaus and Tenenhaus, 2011\].

#'
#' @param X numeric matrix (n by p) which contains the observations on the X variables.
#' @param Y numeric matrix (n by q) which contains the observations on the Y variables.
#'
#' @return A list with the following components of the Canonical Structure:
#' \item{cor }{Canonical correlations.}
#' \item{p_values }{The corresponding p-values for the each of  the canonical correlations.}
#' \item{canonical_weights_X}{The canonical weights for the variables of the dataset X.}
#' \item{canonical_weights_Y}{The canonical weights for the variables of the dataset Y.}
#' \item{canonical_factor_loadings_X}{The inter-set canonical factor loadings for the variables of the dataset X.}
#' \item{canonical_factor_loadings_Y}{The inter-set canonical factor loadings for the variables of the dataset Y.}
#'
#' @references Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis,
#'             Computational Statistics & Data Analysis, Volume 70, 2014, Pages 88-100,
#'             ISSN 0167-9473, https://doi.org/10.1016/j.csda.2013.09.020.
#'
#' Schafer, J; Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation
#' and implications for functional genomics. Statistical Applications in Genetics and Molecular
#' Biology 4:14, Article 32.
#'
#' Tenenhaus, A.; Tenenhaus, M. (2011). Regularized Generalized Canonical Correlation Analysis.
#' Psychometrika 76:2, DOI: 10.1007/S11336-011-9206-8.
#' @author Raul Cruz-Cano
#' @export
#' @examples # Example # 1 Multivariate Normal Data
#' p<-10
#' q<-10
#' n<-50
#' res<-generate_multivariate_normal_sample(p,q,n)
#' X<-res$X
#' Y<-res$Y
#' rownames(X)<-c(1:n)
#' colnames(X)<-c(1:p)
#' colnames(Y)<- c(1:q)
#' my_res<-frcc(X,Y)
#' print(my_res)
#' #Example #2 Soil Specification Data
#' data(soilspec)
#' list_of_units_to_be_used<-sample(1:nrow(soilspec),14)
#' X<- soilspec[list_of_units_to_be_used,2:9]
#' Y<- soilspec[list_of_units_to_be_used,10:15]
#' colnames(X)<-c("H. pubescens", "P. bertolonii", "T. pretense",
#'                "P. sanguisorba", "R. squarrosus", "H. pilosella", "B. media","T. drucei")
#' colnames(Y)<- c("d","P","K","d x P", "d x K","P x K")
#' my_res<-frcc(X,Y)
#' grDevices::dev.new()
#' plot_variables(my_res,1,2)
#' #Example #3 NCI-60 micrRNA Data
#' data("Topoisomerase_II_Inhibitors")
#' data("microRNA")
#' my_res <- frcc(t(microRNA),-1*t(Topoisomerase_II_Inhibitors))
#' for( i in 1:dim(microRNA)[2])
#' {
#'   colnames(microRNA)[i]<-substr(colnames(microRNA)[i], 1, 2)
#' }#end for i
frcc<-function(X,Y)
{
  #Dimensions of the problem
  p<-dim(X)[2]
  #p = #number of variables in X
  q<-dim(Y)[2]
  #q = #number of variables in Y
  n<-dim(X)[1] #or =dim(Y)[1]
  #n = #number of observations
  #not everyting is declared but ut works just fine
  res.frcc<-list(cor=matrix(data=c(0),nrow=1,ncol=q),p_values=matrix(data=c(0),nrow=1,ncol=q),canonical_weights_X=matrix(data=c(0),nrow=p,ncol=q),canonical_weights_Y=matrix(data=c(0),nrow=q,ncol=q))
  #just in case that we receive a datasets that has not been standarized:
  for (i in 1:p)
  { X[,i]<-(X[,i]-mean(X[,i]))/as.numeric(sapply(as.data.frame(X[,i]),stats::sd))
  } #end for i
  for (i in 1:q)
  { Y[,i]<-(Y[,i]-mean(Y[,i]))/as.numeric(sapply(as.data.frame(Y[,i]),stats::sd))
  } #end for i
  #------finding regularized cross-correlation matrix---------------
  S_star_XX <- corpcor::cor.shrink (X, verbose = FALSE) # ADDED VERBOSE = FALSE (AARON)
  S_star_YY <- corpcor::cor.shrink (Y, verbose = FALSE) # ADDED VERBOSE = FALSE (AARON)
  lambda_XY<-off.diagonal.lambda(cbind(X,Y), p,q)
  #shrinking the off-diagonal matrices
  S_star_XY <- stats::cor(X,Y) *(1- lambda_XY)# +lambda_XY*T_XY, but T_XY=null-matrix
  S_star_YX<- t(S_star_XY)
  res.frcc$lambda_XY <- lambda_XY # ADDED LAMBDA_XY TO THE OBJECT (AARON)
  #--------------------calculating the canonical weights    -------------------------------------------------------
  latent_roots<-eigen(solve(S_star_YY)%*%S_star_YX%*%solve(S_star_XX)%*%S_star_XY, only.values = FALSE) #from 1st Eq. page 14 "Understanding..." book
  res.frcc$canonical_weights_Y<- latent_roots$vectors
  #calculating canonical weights for X
  for (i in 1:q)
  {
    res.frcc$canonical_weights_X[,i] <- solve(S_star_XX)%*%S_star_XY %*% res.frcc$canonical_weights_Y[,i]
  }#end for i
  #Assigning names to the variables
  rownames(res.frcc$canonical_weights_X)<-colnames(X)
  rownames(res.frcc$canonical_weights_Y)<-colnames(Y)
  #------Calculating Canonical Correlations Original Method, not as in Clark book (can you believe that the guy from the CCP packae called them canonical coefficients?)---------------------------
  for (i in 1:q)
  {
    U<- as.matrix(X) %*% as.numeric(res.frcc$canonical_weights_X[,i])
    U<-(U-mean(U))/as.numeric(sapply(as.data.frame(U),stats::sd))

    #vertical axis for the ith CCA units
    V<- as.matrix(Y) %*% res.frcc$canonical_weights_Y[,i]
    V<-(V-mean(V))/as.numeric(sapply(as.data.frame(V),stats::sd))

    res.frcc$cor[i]<-stats::cor(U,V)
  }#end for i
  #obtaining p-values
  #------Checking if the correlations are in order, if we have too little data =>they are not
  flag_anomaly_in_correlations<-0
  for(i in 1:(q-1))
  { if (res.frcc$cor[i] < res.frcc$cor[i+1])
  {flag_anomaly_in_correlations <-1}
  }
  if(flag_anomaly_in_correlations > .5)
  {
    warning("Too little data, results might be unreliable")
    res.frcc <- rearrange.frcc(res.frcc)
  }
  #Calculating the p-values
  p_values_struct=0
  #p_values_struct<-CCP::p.asym(res.frcc$cor, n, p, q, tstat = "Wilks") #we can use   p.perm if we are not sure about the normality of the data
  invisible(utils::capture.output(p_values_struct<-CCP::p.asym(res.frcc$cor, n, p, q, tstat = "Wilks"))) # SUPRESS STD. OUTPUT (AARON)
  res.frcc$p_values<- p_values_struct$p.value
  #-----------Calculating Canonical Factor loadings (Statistica) or inter-set correlation coefficients (Gittins)-----------------
  canonical_factor_loadings_X <- matrix(data=c(0),nrow=p,ncol=q)    #page  38 and 39 of Gittins' book  (just not the final formula but the2nd after the definition)
  canonical_factor_loadings_Y <- matrix(data=c(0),nrow=q,ncol=q)
  #S_star_XY<-cor(X,Y)
  #S_star_YX<- t(S_star_XY)
  for (i in 1:q)
  {
    #canonical_factor_loadings_X[,i] <-  S_star_XY %*%  as.matrix(res.frcc$canonical_weights_Y[,i])
    # Our data is already standarized hence our X is z^(x) and Y is z^(y) in Gittins, b= canonical_weights_Y
    #print(dim(as.matrix( res.frcc$canonical_weights_Y[,i] )))
    #print(class(Y))
    #print(dim(as.matrix(Y)))
    canonical_factor_loadings_X[,i] <- stats::cor(as.matrix(X), as.matrix(Y) %*% as.matrix( res.frcc$canonical_weights_Y[,i] ) )
    #canonical_factor_loadings_Y[,i] <-  S_star_YX %*%  as.matrix(res.frcc$canonical_weights_X[,i])
    canonical_factor_loadings_Y[,i] <- stats::cor(as.matrix(Y), as.matrix(X) %*%  as.matrix(res.frcc$canonical_weights_X[,i])  )
  }#end for i
  res.frcc$canonical_factor_loadings_X <- matrix(data=c(0),nrow=p,ncol=q)    #page  38 and 39 of Gittins' book
  res.frcc$canonical_factor_loadings_Y <- matrix(data=c(0),nrow=q,ncol=q)
  #As in [Schafer and Strimmer] top of p.10 we need to take care of value s greater than 1
  for (i in 1:q)
  {
    res.frcc$canonical_factor_loadings_X[,i] <- canonical_factor_loadings_X[,i]
    res.frcc$canonical_factor_loadings_Y[,i] <- canonical_factor_loadings_Y[,i]
  }#end for i
  #Assigning names to the variables
  rownames(res.frcc$canonical_factor_loadings_X)<-colnames(X)
  rownames(res.frcc$canonical_factor_loadings_Y)<-colnames(Y)
  return(res.frcc)
}#end mrc function
