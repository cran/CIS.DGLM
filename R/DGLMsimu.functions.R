#' @title Simulated Data with DGLM model
#'
#' @description This function implements a simulation based on randomly-generated data following a known structure to create a double generalized linear model (DGLM). In particular, it supports the use of interaction terms in the DGLM.
#'
#' @param n.rep The number of repetitions of each environment level. Defaults to 3.
#' @param n.level.env The number of environment variable levels. Defaults to 2.
#' @param n.obs.per.rep A parameter to accommodate repetitions in the data set. Defaults to 150.
#' @param n.loci The number of gene loci in data set. Defaults to 8.
#' @param which.mean.loci A vector that specifies which gene loci are significant in the mean model. Defaults to the vector c(3:4).
#' @param hypo.mean.para A vector that contains the slopes of the gene locations (from which.mean.loci) in the mean model. Defaults to the vector c(1,4).
#' @param incept.mean The intercept for the mean model portion of the DGLM. Defaults to 36.
#' @param which.var.loci A vector that specifies which gene loci are significant in the variance model. Defaults to c(4:5).
#' @param hypo.var.para A vector that contains the slopes of the gene locations (from which.var.loci) in the variance model. Defaults to c(2,5).
#' @param incept.var The intercept for the variance model portion of the DGLM. Defaults to -2.
#' @param which.env.inter.mean Specifies which level of environment significantly interacts in the mean model. This is one of the parameters that supports model-building with interaction terms. Defaults to 2.
#' @param which.loci.inter.mean Specifies which gene location interacts with the environment level in the mean model. This is one of the parameters that supports model-building with interaction terms. Defaults to 7.
#' @param hypo.inter.para.mean The interaction effect between the specified environment level and the specified gene location in the mean model. This is another parameter that supports model-building with interaction terms. Defaults to 2.5.
#' @param which.env.inter.var Specifies which level of environment significantly interacts in the variance model. This is another parameter that supports model-building with interaction terms. Defaults to 2.
#' @param which.loci.inter.var Specifies which gene location interacts with the environment level in the variance model. This is another parameter that supports model-building with interaction terms. Defaults to 8.
#' @param hypo.inter.para.var Interaction effect between the specified environment level and the specified gene location in the variance model. Defaults to 3.8.
#' @param simu.prob.var The probability of each lcoi to be zero or one. Different loci may have different probability and can be adjusted as needed. Defaults to rep(0.5,n.loci).
#' @param rc.val A parameter used to control right-censored data. It sets an upper limit such that if an observation is above a certain range, that observation cannot be included in the data set. Default value is 40.
#' @param ran.seed The random seed to be used. Defaults to NULL.
#'
#' @return  The function returns a data frame containing the data built from the simulation. It provides data for stress values, an environment variable with two levels (0 and 1), and levels for each simulated gene variable.
#' @import dplyr
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' @export
simu.inter.dat.interboth<-function(n.rep=3, n.level.env=2, n.obs.per.rep=150,
                                   n.loci=8, which.mean.loci=c(3:4), hypo.mean.para=c(1,4), incept.mean=36,
                                   which.var.loci=c(4:5),  hypo.var.para=c(2,5), incept.var=-2,
                                   which.env.inter.mean=2,which.loci.inter.mean=7, hypo.inter.para.mean=2.5,
                                   which.env.inter.var =2,which.loci.inter.var =8, hypo.inter.para.var =3.8,
                                   simu.prob.var=rep(0.5, n.loci), rc.val=40,
                                   ran.seed=NULL)
{
  #generate random data with known structure, replicating environment level, and genes.
  # With the provided parameters we can estimate stress value.
  if(is.null(ran.seed)!=1) set.seed(ran.seed);

  n.mean.loci<-ifelse(length(which.mean.loci)==length(hypo.mean.para),length(which.mean.loci),NA)
  n.var.loci<-ifelse(length(which.var.loci)==length(hypo.var.para),length(which.var.loci),NA)

  if (which.env.inter.mean>n.level.env) {print("the which.env.sig should be less than n.level.env"); stop;}
  if (which.env.inter.var >n.level.env) {print("the which.env.sig should be less than n.level.env"); stop;}

  env<-rep(rep(c(1:n.level.env),rep(n.rep,n.level.env)),n.obs.per.rep);
  env.mat<-matrix(0,nrow=n.rep*n.level.env*n.obs.per.rep, ncol=n.level.env);
  for(i in 1:n.level.env)  env.mat[,i]=(env==i);
  var.names.env<-apply(expand.grid('env_var', c(1:n.level.env)), 1, paste, collapse=".");
  colnames(env.mat)<-var.names.env;

  simu.loci<-NULL;
  for(i in 1:n.loci) simu.loci<-cbind(simu.loci,rep(stats::rbinom(n.obs.per.rep, 1,simu.prob.var[i]),rep(n.rep*n.level.env,n.obs.per.rep)))

  var.names.loci<-apply(expand.grid('loci_var', c(1:n.loci)), 1, paste, collapse=".");
  colnames(simu.loci)<-var.names.loci;

  all.mean.para<-rep(0,n.loci);all.mean.para[which.mean.loci]<-hypo.mean.para;
  all.var.para<-rep(0,n.loci);all.var.para[which.var.loci]<-hypo.var.para;


  simu.obs<-cbind(env.mat[,-1],simu.loci); ## the -1 is to get rid of the first level for collinearity.
  colnames(simu.obs)[1:(n.level.env-1)]<-var.names.env[-1];
  ### question: do we need to have similar data structure as the real data set, that is, the same type of genes? 🧬🧬🧬
  ### or, would a gene 🧬🧬🧬 structure(a composition of multiple A' 🅰️ s and B' 🅱️s from gene🧬🧬🧬 loci) be an important thing to consider, or just with gene loci is contributing to the phenotype?
  inter.mat<-NULL;
  for(i in n.level.env:1) for(j in n.loci:1) inter.mat<-cbind(simu.loci[,j]*env.mat[,i],inter.mat); ## stop at 2 is to get rid of the first level for collinearity.

  var.names.inter.mat<-apply(as.matrix(expand.grid(var.names.loci, var.names.env)),1,paste,collapse="*");
  colnames(inter.mat)<-var.names.inter.mat;

  all.inter.para.mean<- all.inter.para.var<-rep(0,n.level.env*n.loci);
  all.inter.para.mean[which.loci.inter.mean+(which.env.inter.mean-1)*n.loci]<-hypo.inter.para.mean;
  all.inter.para.var[which.loci.inter.var +(which.env.inter.var -1)*n.loci]<-hypo.inter.para.var;

  #which.env.inter.mean=3,which.loci.inter.mean=7, hypo.inter.para.mean=2.5,
  inter.pyntp.mean<-inter.mat%*%all.inter.para.mean; # accommodate all the interaction effects
  inter.pyntp.var <-inter.mat%*%all.inter.para.var;

  mean.pyntp<-incept.mean+simu.loci%*%all.mean.para;
  var.pyntp<-sapply(incept.var+simu.loci%*%all.var.para+inter.pyntp.var,FUN=function(x) stats::rnorm(1,sd=sqrt(exp(x))));


  stress<-mean.pyntp+var.pyntp+inter.pyntp.mean # main function to generate the mean model, the variance part and interaction part (both mean and variance)
  simu.obs.df<-data.frame(stress,simu.obs,inter.mat[,-c(1:n.loci)]); #make it a data frame

  return(simu.obs.df);
}
# the above function is to generate the simulated data


# now we analyze the data we generated above by using the forward stepwise selection approach

#' @title Forward Stepwise Selection for Simulated Data
#'
#' @description This function implements the forward stepwise variable selection procedure on the simulated data set generated in simu.inter.dat.interboth. In this function, we utilize a dummy value of "1" when initializing the model to avoid issues with a NULL value when adding variables to the model.
#'
#' @param dat.ana.num12.df A data set filled with data based on a simulation, per the procedure to generate it in the simu.inter.dat.interboth function.
#' @param ouput.name The name of the output file to which the results are to be saved. Defaults to 'out1.txt'.
#' @param num.loop The number of iterations that forward stepwise selection is performed (and hence how many variables will be in the final mean and variance models). Defaults to 10 loops.
#'
#' @return A list with mean and variance mean effects and p-values associated with the coefficients.
#' @import dplyr
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' forward.sel.dglm(test.data)
#' @export
forward.sel.dglm<-function(dat.ana.num12.df, ouput.name='out1.txt', num.loop=10)
{ # dat.ana.num12.df is a dataframe, the response is "stress"
  # num.loop=30 is the max step to have
  sink(ouput.name);
  current.mean.pool.LO<-current.var.pool.LO<-colnames(dat.ana.num12.df)[-1];
  forward.mean.pool<-forward.var.pool<-1; # these are initialized to 1; cannot use NULL, or error is generated
  conflict.set<-vector(mode = "list", length = num.loop);

  for(k in 1:num.loop)
  {
    print(paste('-------------------------------------Running the',k,'th loop-----------------------------------------------' ));
    sig.val.mean<-sig.val.var<-0.05;
    idx.mean<-idx.var<-NULL;
    print(paste('The current mean model is: ', paste('stress ~', paste(forward.mean.pool, collapse="+"),collapse = ''),sep=''))
    for(x in 1:length(current.mean.pool.LO))
    {
      mean.x<-current.mean.pool.LO[x];
      mean.model<- stats::as.formula(paste('stress', paste(c(forward.mean.pool,mean.x), collapse="+"), sep="~"));
      var.model<-switch(is.null(forward.var.pool)+1,stats::as.formula(paste('~', paste(forward.var.pool, collapse="+"), sep="")),stats::as.formula('~ 1'))
      out<-tryCatch( {mod.forward.update<-dglm::dglm( mean.model, dformula=var.model,data=dat.ana.num12.df);
      res.mod.update<-summary(mod.forward.update);
      if (res.mod.update$coefficients%>%last<sig.val.mean) {
        sig.val.mean<-res.mod.update$coefficients%>%last;
        idx.mean<-x;}
      },
      error=function(cond) {print(paste(mean.x, 'causing error when entering the model'));
        return(cond)},
      warning=function(cond) {print(paste(mean.x, 'causing Warning when entering the model'));
        return(cond)}
      );

    }
    if(!is.null(idx.mean)) { print(paste('Adding [',current.mean.pool.LO[idx.mean],'] to mean model, with p-value=',sig.val.mean,sep=''));
      forward.mean.pool<-c(forward.mean.pool,current.mean.pool.LO[idx.mean]);
      current.mean.pool.LO<-current.mean.pool.LO[-idx.mean];
    } ##  else
    print(paste('The current Var model is: ', paste(' ~', paste(forward.var.pool, collapse="+"),collapse = ''),sep=''))

    for(z in 1:length(current.var.pool.LO))
    {
      var.z<-current.var.pool.LO[z];
      mean.model<- stats::as.formula(paste('stress', paste(c(forward.mean.pool), collapse="+"), sep="~"));
      var.model<- stats::as.formula(paste('~', paste(c(forward.var.pool,var.z), collapse="+"), sep=""));
      out<-tryCatch( {mod.forward.update<-dglm::dglm( mean.model, dformula=var.model,data=dat.ana.num12.df);
      res.mod.update<-summary(mod.forward.update);
      if (res.mod.update$dispersion.summary$coefficients%>%last<sig.val.var) {
        sig.val.var<-res.mod.update$dispersion.summary$coefficients%>%last;
        idx.var<-z;}
      },
      error=function(cond) {print(paste(var.z, 'causing error when entering the model'));
        return(cond)},
      warning=function(cond) {print(paste(var.z, 'causing Warning when entering the model'));
        return(cond)}
      );
    }

    if(!is.null(idx.var)) {print(paste('Adding [',current.var.pool.LO[idx.var],'] to var model, with p-value=',sig.val.var,sep=''));
      forward.var.pool<-c(forward.var.pool,current.var.pool.LO[idx.var]);
      current.var.pool.LO<-current.var.pool.LO[-idx.var];
    }

  }

  final.mean.model<- stats::as.formula(paste('stress', paste(c(forward.mean.pool), collapse="+"), sep="~"));
  final.var.model<- stats::as.formula(paste('~', paste(c(forward.var.pool), collapse="+"), sep=""));
  final.mod.forward.update<-dglm::dglm(final.mean.model, dformula=final.var.model,data=dat.ana.num12.df);
  print('-----------------------------Final Model ----------------------------')
  summary(final.mod.forward.update);
  print(final.mean.model);
  print(final.var.model);
  print('--------------------Conflict variable with models ----------------------------')
  print(conflict.set)
  sink(NULL);

  sum.model<-summary(final.mod.forward.update)

  main.coef.eff<-colSums(sum.model$coefficients[,1][-1]*outer(row.names(sum.model$coefficients)[-1],colnames(dat.ana.num12.df)[-1],'=='))
  main.coef.pval<-colSums(sum.model$coefficients[,4][-1]*outer(row.names(sum.model$coefficients)[-1],colnames(dat.ana.num12.df)[-1],'=='))
  main.coef.eff[main.coef.eff==0]<-NA;main.coef.pval[is.na(main.coef.eff)]<-NA;

  var.coef.eff<-colSums(sum.model$dispersion.summary$coefficients[,1][-1]*outer(row.names(sum.model$dispersion.summary$coefficients)[-1],colnames(dat.ana.num12.df)[-1],'=='))
  var.coef.pval<-colSums(sum.model$dispersion.summary$coefficients[,4][-1]*outer(row.names(sum.model$dispersion.summary$coefficients)[-1],colnames(dat.ana.num12.df)[-1],'=='))
  var.coef.eff[var.coef.eff==0]<-NA;var.coef.pval[is.na(var.coef.eff)]<-NA;

  return(list(main.coef.eff,main.coef.pval, var.coef.eff,var.coef.pval))
}

