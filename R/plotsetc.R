#===================draws circles========================================================
#from http://www.r-bloggers.com/circle-packing-with-r/ by Michael Bedward
#' Draws a circle
#'
#' @description Given a center, radius and color, this function draws a circle.
#'
#' @param x X coordinate of the center
#' @param y Y coordinate of the center
#' @param r Radius of the circle
#' @param col Color of the circle
#'
#' @return This function does not return a value, it just draws a circle.
#' @references  http://www.r-bloggers.com/circle-packing-with-r/
#' @keywords Internal
#' @author Michael Bedward
#' @export
custom.draw.circle <- function(x, y, r, col) {
  graphics::lines( cos(seq(0, 2*pi, pi/180)) * r + x, sin(seq(0, 2*pi, pi/180)) * r + y , col=col )
}
#=====================Plot Units=========================================================
#' Plots the experimental units in the Canonical Variates Space
#'
#' @description This function plots the experimental units used in the FRCCA as points in a two-dimensional
#'   plane in which the axis are the canonical variates selected by the user
#'
#' @param X numeric matrix (n by p) which contains the observations on the X variables.
#' @param Y numeric matrix (n by p) which contains the observations on the Y variables.
#' @param res.mrcc List containing a canonical structure provided by the function frcc for the dataset X and Y.
#' @param i Canonical Variate which will be used for the axes (X for horizontal and Y for vertical).
#' @param text_size Character expansion factor for the labels of the experimental units.
#' @param point_size Character expansion factor for the point representing the experimental units.
#'
#' @return This function just creates the units plot. It does not return a value.
#' @references Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis,
#'             Computational Statistics & Data Analysis, Volume 70, 2014, Pages 88-100,
#'             ISSN 0167-9473, https://doi.org/10.1016/j.csda.2013.09.020.
#'
#' @author Raul Cruz-Cano
#' @export
#' @examples
#' #Example: NCI-60 micrRNA Data
#' data("Topoisomerase_II_Inhibitors")
#' data("microRNA")
#' my_res <- frcc(t(microRNA),-1*t(Topoisomerase_II_Inhibitors))
#' for( i in 1:dim(microRNA)[2])
#' {
#'   colnames(microRNA)[i]<-substr(colnames(microRNA)[i], 1, 2)
#' }#end for i
#' grDevices::dev.new()
#' plot_units(t(microRNA),-1*t(Topoisomerase_II_Inhibitors),my_res,1,1,text_size=0.01)
plot_units<-function(X,Y,res.mrcc,i,text_size=.8,point_size=2)
{
  #horizontal axis for the ith CCA units
  U<- as.matrix(X) %*% as.numeric(res.mrcc$canonical_weights_X[,i])
  U<-(U-mean(U))/as.numeric(sapply(as.data.frame(U),stats::sd))
  #vertical axis for the ith CCA units
  V<- as.matrix(Y) %*% res.mrcc$canonical_weights_Y[,i]
  V<-(V-mean(V))/as.numeric(sapply(as.data.frame(V),stats::sd))
  #------------------Dr. Whitmore's graphs------------------------------------------
  #creating names for the axis
  s1<- paste("U",i,sep="")
  s2<-paste("r=",round(res.mrcc$cor[i], digits = 3),sep="")
  s3<-paste("p-value=",round(res.mrcc$p_value[i], digits = 3))
  my_xlab=paste(s1,s2,s3, sep= "       ")
  my_ylab=paste("V",i)
  plot(U,V, ylim=c(min(V)-.5,max(V)+.5),xlim=c(min(U)-.5,max(U)+.5),col="white",xlab=my_xlab, ylab=my_ylab)
  #Plotting   X cex =size, pch= figure  19=circle, col= color
  graphics::points(U,V, cex=point_size,pch=19,col="black")
  calibrate::textxy(U,V, rownames(X) ,cx=text_size,dcol="black")
}
#=============================Plots Variables============================================
#' Plot variables in the Canonical Factor Loadings Space
#'
#' @description This function plots the variables used in the FRCCA as points in a two-dimensional
#' plane in which the axis are the canonical factor loadings selected by the user.
#'
#' @param res.mrcc   List containing a canonical structure provided by the function frcc.
#' @param i Canonical Factor Loadings which will be used as the horizontal axis.
#' @param j  Canonical Factor Loadings  which will be used as the vertical axis.
#' @param inner_circle_radius Radius of the circle which is used to determine
#'       which variables are significant. Only the significant variables will be labeled.
#' @param text_size Character expansion factor for the labels of the variables.
#'
#' @return  This function just creates the variables plot. It does not return a value.
#' @references Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis,
#'             Computational Statistics & Data Analysis, Volume 70, 2014, Pages 88-100,
#'             ISSN 0167-9473, https://doi.org/10.1016/j.csda.2013.09.020.
#'
#' @author Raul Cruz-Cano
#' @export
#'
#' @examples
#' # Example: Multivariate Normal Data
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
#' grDevices::dev.new()
#' plot_variables(my_res,1,2,text_size=1.0)
plot_variables<-function(res.mrcc,i,j,inner_circle_radius=.5,text_size=.8)
{ #Plots
  p <-  nrow(res.mrcc$canonical_factor_loadings_X)
  q <-  nrow(res.mrcc$canonical_factor_loadings_Y)
  first_dimension<-matrix(data=c(0),nrow=p+q,ncol=1)
  first_dimension[1:p,1] <- res.mrcc$canonical_factor_loadings_X[,i]
  first_dimension[(p+1):(p+q),1] <- res.mrcc$canonical_factor_loadings_Y[,i]

  second_dimension<-matrix(data=c(0),nrow=p+q,ncol=1)
  second_dimension[1:p,1] <- res.mrcc$canonical_factor_loadings_X[,j]
  second_dimension[(p+1):(p+q),1] <- res.mrcc$canonical_factor_loadings_Y[,j]
  #creating names for the axis
  my_xlab=paste("CC",i)
  my_ylab=paste("CC",j)
  #no need to find the limits of the plot, they are always 1, because they are correlations
  plot(first_dimension,second_dimension, ylim=c(-1.25,1.25),xlim=c(-1.25,1.25),pch=" ", xlab=my_xlab, ylab=my_ylab)
  #Plotting   X cex =size, pch= figure  17=triangle, col= color ="grey70"
  graphics::points(first_dimension[1:p,],second_dimension[1:p,], cex=2,pch=17,col="grey70")
  #Plotting   X cex =size, pch= figure
  graphics::points(first_dimension[(p+1):(p+q),],second_dimension[(p+1):(p+q),], cex=2,pch=19,col=1)
  #obtining the names of the variates
  points_names<-c(rownames(res.mrcc$canonical_factor_loadings_X),rownames(res.mrcc$canonical_factor_loadings_Y))
  #deleting the names of the non-significan variables
  for (i in 1:(p+q))
  { if(  sqrt((first_dimension[i,1]^2) + (second_dimension[i,1]^2 ) )< inner_circle_radius )
  {
    points_names[i] <- " "
  }
  }#end for i
  #print(points_names)
  calibrate::textxy(first_dimension,second_dimension, points_names ,cx=text_size,dcol="blue")
  custom.draw.circle(0,0,inner_circle_radius,col="black")
  custom.draw.circle(0,0,1.0,col="black")
  #Lines across the Axis
  graphics::segments(-1.45,0 , 1.45,0 )
  graphics::segments(0,-1.45,0 , 1.45 )
}

#==================Generate multivariate Normal Sample==========================================
#' It generates a sample from a multivariate normal distribution function
#'
#' @description It generates a sample from a multivariate normal distribution
#'     function with the cross-covariance matrix described in \[Cruz-Cano et al. 2012\].
#'
#' @param p Number of desired variables in the dataset X.
#' @param q Number of desired variables in the dataset Y.
#' @param n sample size desired.
#'
#' @return A list of n sample units with the values for the variables of the data sets X and Y.
#'
#' @references Cruz-Cano, R.; Lee, M.L.T.; Fast Regularized Canonical Correlation Analysis,
#'             Computational Statistics & Data Analysis, Volume 70, 2014, Pages 88-100,
#'             ISSN 0167-9473, https://doi.org/10.1016/j.csda.2013.09.020.
#'
#' @author Raul Cruz-Cano
#' @export
#'
#' @examples
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
generate_multivariate_normal_sample<-function(p,q,n)
{
  #p = #number of variables in X
  #q = #number of variables in Y
  #n = #number of observations
  #Creating Correlation matrix, this will be the objective of the Sigmas
  if( (p > 6) & (q>6))
  {
    Sigma_Z<-matrix(data=c(0.0),nrow<-(p+q), ncol<-(p+q))   # From Gnz. Journal of Biological Systems17(2):173-199 (2009) page 178
    for (i in 1:(p+q))
    {  for (j in 1:(q+p))
    { if (i==j)  # Filling the diagonal wih 1's
    {
      #Sigma_Z[i,j]<-1.0
      Sigma_Z[i,j]<- 1.0#*sample(random_number,1)
    }
      if( ((i== (p+1)) & (j==1)) || ((j==(p+1)) & (i==1)) ) #establishing that X1 and Y1 have a correlation of .9
      {Sigma_Z[i,j]<- .9
      }
      if( ((i== (p+2)) & (j==2)) || ((j==(p+2)) & (i==2)) ) #establishing that X2 and Y2 have a correlation of .7
      {Sigma_Z[i,j]<- .7 # 1.4
      }
      if( ((i== (p+3)) & (j==3)) || ((j==(p+3)) & (i==3)) ) #establishing that X2 and Y2 have a correlation of .7
      {Sigma_Z[i,j]<- .5 # 1.4
      }
      if( ((i== (p+4)) & (j==4)) || ((j==(p+4)) & (i==4)) ) #establishing that X2 and Y2 have a correlation of .7
      {Sigma_Z[i,j]<- .3 # 1.4
      }
      if( ((i== (p+5)) & (j==5)) || ((j==(p+5)) & (i==5)) ) #establishing that X2 and Y2 have a correlation of .7
      {Sigma_Z[i,j]<- .1 # 1.4
      }
    }#end for j
    } #end for i
    #creating random sample from multivariate normal with covariance matrix created above
    means<- matrix(data=c(0.0),nrow<-1, ncol<-(p+q))   #notice mean =0
    ndata<- MASS::mvrnorm(n,Sigma=Sigma_Z,mu=means)
    X<-ndata[,1:p]
    Y<-ndata[,(p+1):(p+q)]
    res<-list(X=matrix(data=c(0),nrow=n,ncol=p),Y=matrix(data=c(0),nrow=n,ncol=q),Sigma_Z=matrix(data=c(0),nrow=(p+q),ncol=(p+q)))
    res$X<-X
    res$Y<-Y
    res$Sigma_Z<- Sigma_Z
    res
  }#end if
  else{stop("Not enough variables")}
}#end function
