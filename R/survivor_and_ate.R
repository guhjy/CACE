#' Basic Survival Analysis
#'
#' @description Get the effect of a treatment on probability of survival
#' and on the outcome variable among the survivors. For the first, we
#' take the ATE function from npcausal (Edward Kennedy's package). For
#' the second, we adapt this ATE code slightly.
#'
#' @param y is the outcome vector
#' @param a is the treatment vector (dichotomous)
#' @param t is an indicator for survival
#' @param x is a dataframe of covariates
#'
#' @return Two lists, each containing the following components:
#' \item{res}{ estimates/SEs/CIs/p-values for population means and relevant contrasts.}
#' \item{nuis}{ subject-specific estimates of nuisance functions (i.e., propensity score and outcome regression) }
#' \item{ifvals}{ matrix of estimated influence function values.}
#'
#' @example
#' n <- 1000; x <- matrix(rnorm(n*5),nrow=n)
#' a <- rbinom(n, 1,.3); t <- rbinom(n,1,.8)
#' y <- rnorm(n, 1+.5*a, 1); y[t==0] <- 0
#' res <- survival(y,a,t,x)

survival <- function(y,a,t,x,nsplits=2,
                 sl.lib=c("SL.earth","SL.gam","SL.glm","SL.glm.interaction","SL.mean","SL.rpart")){

eff.on.survival <- ate(t,a,x, nsplits = nsplits, sl.lib = sl.lib)
eff.on.survivors <- surv2(y,a,t,x, nsplits = nsplits, sl.lib = sl.lib)

return(list(eff.on.survival = eff.on.survival,
            eff.on.survivors = eff.on.survivors))

}



