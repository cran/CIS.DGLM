#' @title Forward Stepwise Selection for Real Data
#'
#' @description This function implements the forward stepwise variable selection procedure on a real data set. It utilizes the dglm function from the dglm packages to build the model and helps to account for more complex situations such as convergence issues with dglm and interaction terms in the model. In this function, we utilize a dummy value of "1" when initializing the model to avoid issues with a NULL value when adding variables to the model.
#'
#' @param dat.ana.num12.df The data set to be used to build the DGLM.
#' @param ouput.name The name of the output file to which the results will be saved.
#' @param num.loop The number of iterations that forward stepwise selection is performed (and hence how many variables will be in the final mean and variance models). Defaults to 10 iterations.
#' @param typ.err  Type 1 error. The default value is 0.05.
#'
#' @return A data frame with mean and variance mean effects and p-values associated with the coefficients for each loop. The function also produces a text file containing the model-building information at each stage of the loop (i.e. variables causing errors or warnings, the state of the model at each iteration, etc.).
#'
#' @examples
#' library(dplyr)
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' forward.sel.dglm.real(test.data)
#' unlink(c('out1.txt'))
#' @export
forward.sel.dglm.real<-function(dat.ana.num12.df, ouput.name='out1.txt', num.loop=10,typ.err=0.05)
{
  sink(ouput.name);
  current.mean.pool.LO<-current.var.pool.LO<-colnames(dat.ana.num12.df)[-1];
  forward.mean.pool<-forward.var.pool<-1;  # these are initialized to 1; cannot use NULL, or error is generated
  conflict.set<-vector(mode = "list", length = num.loop);
  for(k in 1:num.loop)
  {
    print(paste('-------------------------------------Running the',k,'th loop-----------------------------------------------' ));
    sig.val.mean<-sig.val.var<-typ.err;
    idx.mean<-idx.var<-NULL;
    print(paste('The current mean model is: ', paste('stress ~', paste(forward.mean.pool, collapse="+"),collapse = ''),sep=''))
    for(x in 1:length(current.mean.pool.LO))
    {
      #print(forward.mean.pool)
      mean.x<-current.mean.pool.LO[x];
      mean.model<- stats::as.formula(paste('stress ', paste(c(forward.mean.pool,mean.x), collapse="+"), sep="~"));
      var.model<-switch(is.null(forward.var.pool)+1,stats::as.formula(paste('~', paste(forward.var.pool, collapse="+"), sep="")),stats::as.formula('~ 1'))
      out<-tryCatch( {mod.forward.update<-dglm::dglm( mean.model, dformula=var.model,data=dat.ana.num12.df);
      res.mod.update<-summary(mod.forward.update);
      if (res.mod.update$coefficients%>%last<sig.val.mean) {
        sig.val.mean<-res.mod.update$coefficients%>%last;
        idx.mean<-x;}
      },
      error=function(cond) {print(paste(mean.x, 'causing error when entering the model'));
        return(cond)},
      warning=function(cond) {print(paste(mean.x, 'causing error when entering the model'));
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
      mean.model<- stats::as.formula(paste('stress ', paste(forward.mean.pool, collapse="+"), sep="~"));
      var.model<- stats::as.formula(paste('~', paste(c(forward.var.pool,var.z), collapse="+"), sep=""));
      out<-tryCatch( {mod.forward.update<-dglm::dglm( mean.model, dformula=var.model,data=dat.ana.num12.df);
      res.mod.update<-summary(mod.forward.update);
      if (res.mod.update$dispersion.summary$coefficients%>%last<sig.val.var) {
        sig.val.var<-res.mod.update$dispersion.summary$coefficients%>%last;
        idx.var<-z;}
      },
      error=function(cond) {print(paste(var.z, 'causing error when entering the model'));
        return(cond)},
      warning=function(cond) {print(paste(var.z, 'causing error when entering the model'));
        return(cond)}
      );
    }

    if(!is.null(idx.var)) {print(paste('Adding [',current.var.pool.LO[idx.var],'] to var model, with p-value=',sig.val.var,sep=''));
      forward.var.pool<-c(forward.var.pool,current.var.pool.LO[idx.var]);
      current.var.pool.LO<-current.var.pool.LO[-idx.var];
    }


  }

  final.mean.model<- stats::as.formula(paste('stress ', paste(forward.mean.pool, collapse="+"), sep="~"));
  final.var.model<- stats::as.formula(paste('~', paste(forward.var.pool, collapse="+"), sep=""));
  final.mod.forward.update<-dglm::dglm(final.mean.model, dformula=final.var.model,data=dat.ana.num12.df);
  print('-----------------------------Final Model ----------------------------')
  summary(final.mod.forward.update);
  print(final.mean.model);
  print(final.var.model);
  print('--------------------Conflict variable with models ----------------------------')
  print(conflict.set)
  sink(NULL);

  sum.model<-summary(final.mod.forward.update)

  main.coef.eff<-colSums(sum.model$coefficients[,1][-1]*outer(forward.mean.pool[-1],colnames(dat.ana.num12.df)[-1],'=='))
  main.coef.pval<-colSums(sum.model$coefficients[,4][-1]*outer(forward.mean.pool[-1],colnames(dat.ana.num12.df)[-1],'=='))
  main.coef.eff[main.coef.eff==0]<-NA;main.coef.pval[is.na(main.coef.eff)]<-NA;

  var.coef.eff<-colSums(sum.model$dispersion.summary$coefficients[,1][-1]*outer(forward.var.pool[-1],colnames(dat.ana.num12.df)[-1],'=='))
  var.coef.pval<-colSums(sum.model$dispersion.summary$coefficients[,4][-1]*outer(forward.var.pool[-1],colnames(dat.ana.num12.df)[-1],'=='))
  var.coef.eff[var.coef.eff==0]<-NA;var.coef.pval[is.na(var.coef.eff)]<-NA;
  res<-rbind(main.coef.eff,main.coef.pval, var.coef.eff,var.coef.pval)%>%data.frame;
  colnames(res)<-colnames(dat.ana.num12.df)[-1];
  rownames(res)<-c('main.eff','main.pval','var.eff','var.pval')
  return(res)
}
### the function above performs the forward selection and produces a text file with the significant variables for mean and variance models and their corresponding pvalues
