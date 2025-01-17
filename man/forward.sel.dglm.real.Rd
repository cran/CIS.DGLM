% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forselhyb.functions.R
\name{forward.sel.dglm.real}
\alias{forward.sel.dglm.real}
\title{Forward Stepwise Selection for Real Data}
\usage{
forward.sel.dglm.real(
  dat.ana.num12.df,
  ouput.name = "out1.txt",
  num.loop = 10,
  typ.err = 0.05
)
}
\arguments{
\item{dat.ana.num12.df}{The data set to be used to build the DGLM.}

\item{ouput.name}{The name of the output file to which the results will be saved.}

\item{num.loop}{The number of iterations that forward stepwise selection is performed (and hence how many variables will be in the final mean and variance models). Defaults to 10 iterations.}

\item{typ.err}{Type 1 error. The default value is 0.05.}
}
\value{
A data frame with mean and variance mean effects and p-values associated with the coefficients for each loop. The function also produces a text file containing the model-building information at each stage of the loop (i.e. variables causing errors or warnings, the state of the model at each iteration, etc.).
}
\description{
This function implements the forward stepwise variable selection procedure on a real data set. It utilizes the dglm function from the dglm packages to build the model and helps to account for more complex situations such as convergence issues with dglm and interaction terms in the model. In this function, we utilize a dummy value of "1" when initializing the model to avoid issues with a NULL value when adding variables to the model.
}
\examples{
library(dplyr)
test.data <- simu.inter.dat.interboth(n.rep = 3, n.obs.per.rep = 15, ran.seed = 1)
forward.sel.dglm.real(test.data)
unlink(c('out1.txt'))
}
