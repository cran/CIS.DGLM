#' @title Mean Stress
#'
#' @description This function provides the mean stress among As and Bs, corresponding to different environment levels, for a list of variables.
#'
#' @param dataset Data set to be utilized. For the purposes of this function, the binary values in each variable are considered to be -1 and 1.
#' @param var A list of variables. If using variables from a DGLM, use the variables from the mean model (or, if trying to find intervals for use in the plotting functions draw.crossplots and draw.squareplots, use all variables in both mean and variance models).
#' @param stress_variable Name of the variable with the stress values.
#' @param output.name Name of the output file to which to save the outputs. Defaults to 'mean_stress.txt'.
#' @param use.output A binary variable to indicate whether the output is automatically saved to an external text file. Defaults to TRUE. If FALSE, the output will not be saved to a file.
#' @param bin.levels A list that provides the binary values utilized in the dataset. Defaults to c(0,1), indicating that 0 and 1 are used as the binary outcomes; can also be 1, -1. List the value for the "A" environment level first, then the value for the "B" environment level.
#' @return Produces a data frame with three columns: var, AvgB, and AvgA. These provide the variable and its corresponding mean stress values for As and Bs, corresponding to different environment levels.
#'
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 2, n.obs.per.rep = 15, ran.seed = 1)
#' variables <- colnames(test.data[-1])
#' mean_stress(test.data, variables, 'stress', use.output = FALSE)
#' @export
mean_stress<-function(dataset,var,stress_variable,output.name='mean_stress.txt',use.output=TRUE, bin.levels=c(0,1)){
  if (use.output){
    sink(output.name);
  }
  out = data.frame(var=character(), AvgB=double(), AvgA=double())
  for (i in var){
    v.txt<-paste('dataset',i,sep='$')
    v<-eval(parse(text = v.txt));
    v2.txt<-paste('dataset',stress_variable,sep='$')
    v2<-eval(parse(text = v2.txt));
    d = data.frame(v,v2)
    AvgB = mean((d[d[,'v']==bin.levels[2],'v2'])) # avg stress among Bs (or 1s)
    AvgA = mean((d[d[,'v']==bin.levels[1],'v2'])) # avg stress among As (or -1s)
    out[nrow(out)+1,'var'] <- i
    out[nrow(out),'AvgB'] <- AvgB
    out[nrow(out),'AvgA'] <- AvgA
  }
  #if (use.output){
  #  sink();
  #}
  return(out)
}
### the function above takes the dataset, a list of variables, and the stress variable and produces a data frame with the
###     average stress of As and average stress of Bs for each variable listed in var


#' @title Standard Deviation Stress
#'
#' @description This function provides the mean stress among As and Bs, corresponding to different environment levels, for a list of variables.
#'
#' @param dataset Data set to be utilized. For the purposes of this function, the binary values in each variable are considered to be -1 and 1.
#' @param var A list of variables. If using variables from a DGLM, use the variables from the variance model (or, if trying to find intervals for use in the plotting functions draw.crossplots and draw.squareplots, use all variables in both mean and variance models).
#' @param stress_variable Name of the variable with the stress values.
#' @param output.name Name of the output file to which to save the outputs. Defaults to 'sd_stress.txt'.
#' @param use.output A binary variable to indicate whether the output is automatically saved to an external text file. Defaults to TRUE. If FALSE, the output will not be saved to a file.
#' @param bin.levels A list that provides the binary values utilized in the dataset. Defaults to c(0,1), indicating that 0 and 1 are used as the binary outcomes; can also be 1, -1. List the value for the "A" environment level first, then the value for the "B" environment level.
#' @return Produces a data frame with three columns: var, sd1, and sdneg1. These provide the variable and its corresponding standard deviation of stress values for As and Bs, corresponding to different environment levels.
#'
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' variables <- colnames(test.data[-1])
#' sd.stress(test.data, variables, 'stress', use.output = FALSE)
#' @export
sd.stress<-function(dataset,var,stress_variable,output.name='sd_stress.txt',use.output=TRUE, bin.levels=c(0,1)){
  if (use.output){
    sink(output.name);
  }
  out = data.frame(var=character(), sdB=double(), sdA=double())
  for (i in var){
    v.txt<-paste('dataset',i,sep='$')
    v<-eval(parse(text = v.txt));
    v2.txt<-paste('dataset',stress_variable,sep='$')
    v2<-eval(parse(text = v2.txt));
    d = data.frame(v,v2)
    sdB = stats::sd((d[d[,'v']==bin.levels[2],'v2'])) # sd for stress among Bs (or 1s)
    sdA = stats::sd((d[d[,'v']==bin.levels[1],'v2'])) # sd for stress among As (or -1s)
    out[nrow(out)+1,'var'] <- i
    out[nrow(out),'sdB'] <- sdB
    out[nrow(out),'sdA'] <- sdA
  }
  #if (use.output){
  #  sink();
  #}
  return(out)
}
### the function above takes the dataset, a list of variables, and the stress variable and produces a dataframe with the
###     standard deviation of stress of As and standard deviation of stress of Bs for each variable listed in var


#' @title Bootstrap
#'
#' @description This function implements a custom bootstrapping procedure that utilizes bootstrapping to estimate mean and SD of stress between two environment states (A and B).
#'
#' @param dataset Data set to be utilized.
#' @param n.boot Number of bootstraps to perform. Defaults to 10^5.
#' @param variables List of variables from mean and variance models in DGLM.
#' @param stress_variable Name of the variable with the stress values.
#' @param alpha Significance level by which to determine the confidence intervals for the bootstrap estimates. Defaults to 0.05, thus creating the 95 percent confidence intervals.
#' @param ran.seed Random seed value for generating different random bootstrap samples.]
#' @return Lists with confidence intervals for the bootstrap estimations for average stress in As and Bs of variables in mean model and confidence intervals for the bootstrap estimations of standard deviation of stress in As and Bs of variables in variance model.
#'
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' variables <- colnames(test.data[-1])
#' bootstrap(test.data, n.boot=100,variables, 'stress')
#' unlink(c('bootstrap mean A stress.txt','bootstrap mean B stress.txt',
#' 'bootstrap sd A stress.txt', 'bootstrap sd B stress.txt'))
#' @export
bootstrap <- function(dataset, n.boot=10^5, variables, stress_variable, alpha=0.05, ran.seed=12345){
  n<-dim(dataset)[1];
  set.seed(ran.seed);
  boot.mean.B<-boot.mean.A<-boot.var.A<-boot.var.B<-NULL;
  i<-1;
  out.loop1<-vector(mode = "list", length = n.boot)
  out.loop2<-vector(mode = "list", length = n.boot)
  while(i<=n.boot)
  {
    boot.index<-sample(1:n,replace=TRUE)
    data.idx <- dataset[boot.index,]
    out.loop1[[i]]<-mean_stress(dataset=data.idx,variables,stress_variable,use.output=FALSE);
    out.loop2[[i]]<-sd.stress(dataset=data.idx,variables,stress_variable,use.output=FALSE);
    # conditional check to ensure that there are no NA's in the final output
    if (!any(is.na(out.loop1[[i]][2])) && !any(is.na(out.loop1[[i]][3])) && !any(is.na(out.loop2[[i]][2])) && !any(is.na(out.loop1[[i]][3]))) {
      # means, separated by A's and B's
      boot.mean.B<-rbind(boot.mean.B,out.loop1[[i]][[2]]);
      boot.mean.A<-rbind(boot.mean.A,out.loop1[[i]][[3]]);
      # SD's, separated by A's and B's
      boot.var.B<- rbind(boot.var.B, out.loop2[[i]][[2]]);
      boot.var.A<- rbind(boot.var.A, out.loop2[[i]][[3]]);
      i<-i+1;
    }
  }
  colnames(boot.mean.A)<- colnames(boot.mean.B) <- colnames(boot.var.A) <- colnames(boot.var.B) <- variables
  mean_stressA <- apply(boot.mean.A,2,stats::quantile,c(alpha/2,1-(alpha/2)))
  mean_stressB <- apply(boot.mean.B,2,stats::quantile,c(alpha/2,1-(alpha/2)))
  var_stressA <- apply(boot.var.A,2,stats::quantile,c(alpha/2,1-(alpha/2)))
  var_stressB <- apply(boot.var.B,2,stats::quantile,c(alpha/2,1-(alpha/2)))

  utils::write.table(mean_stressA, 'bootstrap mean A stress.txt', row.names = T, col.names = T)
  utils::write.table(mean_stressB, 'bootstrap mean B stress.txt', row.names = T, col.names = T)
  utils::write.table(var_stressA, 'bootstrap sd A stress.txt', row.names = T, col.names = T)
  utils::write.table(var_stressB, 'bootstrap sd B stress.txt', row.names = T, col.names = T)

  return (list(mean_stressA, mean_stressB, var_stressA, var_stressB))
}
