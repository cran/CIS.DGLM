#' @title Draw Crossplots
#'
#' @description This function draws crossplots for As and Bs in each variable in the mean and variance models with the Mean Estimate vs Standard Deviation Estimate.
#'
#' @param fn.mean.A Enter file name of file with confidence intervals of mean stress, environment level A data. Can be hybrid, inbred, or full data set. This file needs to be obtained from the bootstrap function, run with the desired data set.
#' @param fn.mean.B Enter file name of file with confidence intervals of mean stress, environment level B data. Can be hybrid, inbred, or full data set. This file needs to be obtained from the bootstrap function, run with the desired data set.
#' @param fn.sd.A Enter file name of file with confidence intervals of SD stress, environment level A data. Can be hybrid, inbred, or full data set. This file needs to be obtained from the bootstrap function, run with the desired data set.
#' @param fn.sd.B Enter file name of file with confidence intervals of SD stress, environment level B data. Can be hybrid, inbred, or full data set. This file needs to be obtained from the bootstrap function, run with the desired data set.
#' @param fn.pe.mean Enter file name of file with point estimates of mean for each gene (both A and B environment levels present). Can be hybrid, inbred, or full data set. This file needs to be obtained from the mean_stress function, run with the desired data set.
#' @param fn.pe.sd Enter file name of file with point estimates of SD for each gene (both A and B environment levels present). Can be hybrid, inbred, or full data set. This file needs to be obtained from the sd.stress function, run with the desired data set.
#' @param variables List of variables from mean and variance models. Mean variables need to be listed first, then variance variables.
#' @param ishybrid Indicates the type of the data set being examined. You can use 'Hybrid', 'Inbred', "All", etc.
#' @param num.vars Number of variables per model. Used to ascertain if a variable falls in the mean or the variance model.
#' @return There is no return for this function; it prints crossplots for each of the variables listed in the parameter 'variables.'
#' @importFrom utils read.table
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' variables <- colnames(test.data[-1])
#' mean_stress(test.data, variables, 'stress')
#' sink();
#' sd.stress(test.data, variables, 'stress')
#' sink();
#' plot_vars <- c("loci_var.4","loci_var.7.env_var.2","loci_var.3",
#' "loci_var.5","loci_var.8.env_var.2","loci_var.4")
#' bootstrap(test.data, n.boot=100,variables, 'stress')
#' draw.crossplots('bootstrap mean A stress.txt','bootstrap mean B stress.txt',
#' 'bootstrap sd A stress.txt', 'bootstrap sd B stress.txt', 'mean_stress.txt',
#' 'sd_stress.txt', plot_vars, 'All',3)
#' unlink(c('bootstrap mean A stress.txt','bootstrap mean B stress.txt',
#' 'bootstrap sd A stress.txt', 'bootstrap sd B stress.txt',
#' 'mean_stress.txt', 'sd_stress.txt'))
#' @export
draw.crossplots <- function(fn.mean.A, fn.mean.B, fn.sd.A, fn.sd.B, fn.pe.mean, fn.pe.sd, variables, ishybrid, num.vars) {
  hyb.mean.A <- utils::read.table(fn.mean.A, header = T)
  hyb.mean.B <- utils::read.table(fn.mean.B, header = T)
  hyb.sd.A <- utils::read.table(fn.sd.A, header = T)
  hyb.sd.B <- utils::read.table(fn.sd.B, header = T)

  # data frame with mean and sd point estimates
  # note: -1 = A; 1 = B
  # to access A & B in the data table, you need a table in the format as that produced by mean_stress and sd.stress in job2.diff.stress
  # A is -1 & is in the third column; B is 1 & is in the second column
  pe.mean <- read.table(fn.pe.mean, header = T)
  pe.sd <- read.table(fn.pe.sd, header = T)

  mean_estimate <- sd_estimate <- xmin <- xmax <- ymin <- ymax <- NULL
  for (i in seq(1, length(variables), by=1)){
    df <- data.frame(mean_estimate = c(pe.mean[i,3],pe.mean[i,2]), sd_estimate = c(pe.sd[i,3],pe.sd[i,2]),
                     xmin = c(hyb.mean.A[1,i],hyb.mean.B[1,i]), xmax = c(hyb.mean.A[2,i], hyb.mean.B[2,i]),
                     ymin = c(hyb.sd.A[1,i],hyb.sd.B[1,i]), ymax = c(hyb.sd.A[2,i],hyb.sd.B[2,i]))

    # create crossplot
    # different options based on whether gene is in mean model or variance model
    if (i <= num.vars) {
      print(ggplot2::ggplot(data = df,ggplot2::aes(x = mean_estimate,y = sd_estimate, colour = c('A', 'B'))) +
              ggplot2::geom_point(size=3) +
              ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2") + ggplot2::theme_bw() +
              ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin,ymax = ymax),width=.1,size=1.2) +
              ggplot2::geom_errorbarh(ggplot2::aes(xmin = xmin,xmax = xmax),height=.1,size=1.2) +
              ggplot2::ggtitle(paste("Stress Mean and SD by Env Level,", variables[i], '(mean)', ishybrid, "Data")) + ggplot2::xlab("Mean Estimate, 95% CI") + ggplot2::ylab("SD Estimate, 95% CI") + ggplot2::labs(color=variables[i]))
    }
    else if (i > num.vars) {
      print(ggplot2::ggplot(data = df,ggplot2::aes(x = mean_estimate,y = sd_estimate, colour = c('A', 'B'))) +
              ggplot2::geom_point(size=3) +
              ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2") + ggplot2::theme_bw() +
              ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin,ymax = ymax),width=.1,size=1.2) +
              ggplot2::geom_errorbarh(ggplot2::aes(xmin = xmin,xmax = xmax),height=.1,size=1.2) +
              ggplot2::ggtitle(paste("Stress Mean and SD by Env Level,", variables[i], '(var)', ishybrid, "Data")) + ggplot2::xlab("Mean Estimate, 95% CI") + ggplot2::ylab("SD Estimate, 95% CI") + ggplot2::labs(color=variables[i]))
    }
  }
}




#' @title Draw Square Plots
#'
#' @description This function draws square plots for As and Bs in each variable in the mean and variance models with the Mean Estimate vs Standard Deviation Estimate
#'
#' @param fn.mean.A file name of file with confidence intervals of mean stress, environment level A data. Can be either hybrid or inbred data.
#' @param fn.mean.B file name of file with confidence intervals of mean stress, environment level B data. Can be either hybrid or inbred data.
#' @param fn.sd.A file name of file with confidence intervals of SD stress, environment level A data. Can be either hybrid or inbred data.
#' @param fn.sd.B file name of file with confidence intervals of SD stress, environment level B data. Can be either hybrid or inbred data.
#' @param fn.pe.mean file name of file with point estimates of mean for each gene (both A and B environment levels present). Can be either hybrid or inbred data.
#' @param fn.pe.sd file name of file with point estimates of SD for each gene (both A and B environment levels present). Can be either hybrid or inbred data.
#' @param variables list of variables from mean and variance models. Mean vars needs to be listed first, then variance vars.
#' @param ishybrid indicates the type of the data set being examined. Choose 'Hybrid' or 'Inbred' or "All".
#' @param num.vars number of variables per model. Used to ascertain if a variable falls in the mean or the variance model.
#' @return There is no return for this function; it prints square plots for each of the variables listed in the parameter 'variables'.
#' @examples
#' test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
#' variables <- colnames(test.data[-1])
#' mean_stress(test.data, variables, 'stress')
#' sink();
#' sd.stress(test.data, variables, 'stress')
#' sink();
#' bootstrap(test.data, n.boot=100,variables, 'stress')
#' plot_vars <- c("loci_var.4","loci_var.7.env_var.2","loci_var.3",
#' "loci_var.5","loci_var.8.env_var.2","loci_var.4")
#' draw.squareplots('bootstrap mean A stress.txt','bootstrap mean B stress.txt',
#' 'bootstrap sd A stress.txt', 'bootstrap sd B stress.txt', 'mean_stress.txt',
#' 'sd_stress.txt', plot_vars, 'All', 3)
#' unlink(c('bootstrap mean A stress.txt','bootstrap mean B stress.txt',
#' 'bootstrap sd A stress.txt', 'bootstrap sd B stress.txt',
#' 'mean_stress.txt', 'sd_stress.txt'))
#' @export
draw.squareplots <- function(fn.mean.A, fn.mean.B, fn.sd.A, fn.sd.B, fn.pe.mean, fn.pe.sd, variables, ishybrid, num.vars) {
  hyb.mean.A <- utils::read.table(fn.mean.A, header = T)
  hyb.mean.B <- utils::read.table(fn.mean.B, header = T)
  hyb.sd.A <- utils::read.table(fn.sd.A, header = T)
  hyb.sd.B <- utils::read.table(fn.sd.B, header = T)

  # data frame with mean and sd point estimates
  # note: -1 = A; 1 = B
  # to access A & B in the data table, you need a table in the format as that produced by mean_stress and sd.stress in job2.diff.stress
  # A is -1 & is in the third column; B is 1 & is in the second column
  pe.mean <- read.table(fn.pe.mean, header = T)
  pe.sd <- read.table(fn.pe.sd, header = T)

  mean_estimate <- sd_estimate <- xmin <- xmax <- ymin <- ymax <- NULL
  for (i in seq(1, length(variables), by=1)){
    df <- data.frame(mean_estimate = c(pe.mean[i,3],pe.mean[i,2]), sd_estimate = c(pe.sd[i,3],pe.sd[i,2]),
                     xmin = c(hyb.mean.A[1,i],hyb.mean.B[1,i]), xmax = c(hyb.mean.A[2,i], hyb.mean.B[2,i]),
                     ymin = c(hyb.sd.A[1,i],hyb.sd.B[1,i]), ymax = c(hyb.sd.A[2,i],hyb.sd.B[2,i]))

    # create squareplots
    # different options based on whether gene is in mean model or variance model
    if (i <= num.vars) {
      print(ggplot2::ggplot(data = df,ggplot2::aes(x = mean_estimate,y = sd_estimate, colour = c('A', 'B'))) +
              ggplot2::geom_point(size=3) +
              ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2") + ggplot2::theme_bw() +
              ggplot2::geom_rect(ggplot2::aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax),alpha=0,color=c('chartreuse','lightcoral'),size=1) +
              ggplot2::ggtitle(paste("Stress Mean and SD by Env Level,", variables[i], '(mean)', ishybrid, "Data")) + ggplot2::xlab("Mean Estimate, 95% CI") + ggplot2::ylab("SD Estimate, 95% CI") + ggplot2::labs(color=variables[i]))
    }
    else if (i > num.vars) {
      print(ggplot2::ggplot(data = df,ggplot2::aes(x = mean_estimate,y = sd_estimate, colour = c('A', 'B'))) +
              ggplot2::geom_point(size=3) +
              ggplot2::scale_colour_brewer(type = "qual", palette = "Dark2") + ggplot2::theme_bw() +
              ggplot2::geom_rect(ggplot2::aes(xmin = xmin,xmax = xmax,ymin = ymin,ymax = ymax),alpha=0,color=c('chartreuse','lightcoral'),size=1) +
              ggplot2::ggtitle(paste("Stress Mean and SD by Env Level,", variables[i], '(var)', ishybrid, "Data")) + ggplot2::xlab("Mean Estimate, 95% CI") + ggplot2::ylab("SD Estimate, 95% CI") + ggplot2::labs(color=variables[i]))
    }
  }
}
