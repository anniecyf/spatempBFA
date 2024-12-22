###Function to get model fit diagnostics given a FixedLbfa or varyLjBFA object
#'
#' diagnostics
#'
#' Calculates diagnostic metrics based on outputs from a \code{\link{FixedLbfa}} or \code{\link{varyLjBFA}} model.
#'
#' @param object A \code{\link{FixedLbfa}} or \code{\link{varyLjBFA}} or \code{\link{FixedLbfaVAR1}} or \code{\link{VAR1varyLjBFA}}
#'  model object for which diagnostics are desired from.
#'
#' @param diags A vector of character strings indicating which diagnostics to compute.
#'  Options include: Deviance Information Criterion ("dic"), d-infinity ("dinf"), 
#'  Posterior Mean Information Criteria ("meanIC"), and Watanabe-Akaike Information Criterion ("waic"). 
#'  At least one option must be included.
#'  Note: The probit model cannot compute the DIC or WAIC diagnostics due to computational
#'  issues regarding the multivariate normal CDF.
#'
#' @param keepDeviance A logical indicating whether the posterior deviance distribution is returned (default = FALSE).
#'
#' @param keepPPD A logical indicating whether the posterior predictive distribution
#'  at each observed location is returned (default = FALSE).
#'
#' @param Verbose A boolean logical indicating whether progress should be output (default = TRUE).
#'
#' @param seed An integer value used to set a seed for the random number generator (default = 54).
#'  
#' @details To assess model fit, DIC, d-infinity, meanIC, and WAIC are utilized. DIC is based on the
#'  deviance statistic and penalizes the complexity of a model with an effective
#'  number of parameters estimate pD (Spiegelhalter et al 2002). The d-infinity posterior
#'  predictive measure is an alternative diagnostic tool, where d-infinity=P+G.
#'  The G term decreases as goodness of fit increases, and P, the penalty term, inflates
#'  as the model becomes over-fit. Hence, small values of both of these terms and thus small
#'  values of d-infinity are desirable (Gelfand and Ghosh 1998). The meanIC metric consists of 2 measures 
#'  (one for precision and one for accuracy)
#'  constructed from the posterior predicted response variables for all (i,o,t) at all MCMC iterations.
#'  WAIC is invariant to parametrization and is asymptotically equal to Bayesian cross-validation
#'  (Watanabe 2010). WAIC = -2 * (lppd - p_waic_2). Where lppd is the log pointwise
#'  predictive density and p_waic_2 is the estimated effective number of parameters
#'  based on the variance estimator from Vehtari et al. 2016. (p_waic_1 is the mean
#'  estimator).
#'
#' @return \code{diagnostics} returns a list containing the diagnostics requested and
#'  possibly the deviances corresponding to all kept MCMC iterations.
#'
#' @references Yifan Cheng
#' @references Berchuck, S. I., Janko, M., Medeiros, F. A., Pan, W., & Mukherjee, S. (2021). Bayesian Non-Parametric Factor Analysis for Longitudinal Spatial Surfaces. Bayesian Analysis, 17(2), 1–30.
#' @references Gelfand, A. E., & Ghosh, S. K. (1998). Model choice: a minimum posterior predictive loss approach. Biometrika, 1-11.
#' @references Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 64(4), 583-639.
#' @references Vehtari, A., Gelman, A., & Gabry, J. (2016). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing, 1-20.
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11(Dec), 3571-3594.
#'
#' @export
diagnostics <- function(object, diags = c("dic", "dinf", "meanIC", "waic"), keepDeviance = FALSE, keepPPD = FALSE, 
                        Verbose = TRUE, seed = 29) {
  if (missing(object)) stop('"object" is missing')
  if ( ! ( is.FixedLbfa(object) | is.varyLjBFA(object) | is.FixedLbfaVAR1(object) | is.VAR1varyLjBFA(object) ) ) {
    stop('"object" must be of class FixedLbfa or varyLjBFA or FixedLbfaVAR1 or VAR1varyLjBFA')
  }
  if ( !all(diags %in% c("dic", "dinf", "meanIC", "waic")) | is.null(diags) ) stop('the vector "diags" must be a non-empty subset of c("dic", "dinf", "meanIC", "waic")')
  if (!is.logical(keepDeviance)) stop('"keepDeviance" must be a logical')
  if (!is.logical(keepPPD)) stop('"keepPPD" must be a logical')
  if (!is.logical(Verbose)) stop('"Verbose" must be a logical')
  DatObj <- object$datobj
  set.seed(seed)
  NKeep <- dim(object$lambda)[1]
  M <- DatObj$M
  O <- DatObj$O
  C <- DatObj$C
  K <- DatObj$K
  Nu <- DatObj$Nu
  EyeNu <- diag(Nu)
  FamilyInd <- DatObj$FamilyInd
  N <- M * O * Nu
  YObserved <- DatObj$YStar
  X <- DatObj$X
  
  ###Construct the parameter object
  Para <- list()
  Para$Lambda <- object$lambda
  Para$Eta <- object$eta
  if (DatObj$P > 0) Para$Beta <- object$beta else Para$Beta <- matrix(0, ncol = 0, nrow = NKeep)
  if (is.null(object$sigma2)) Para$Sigma2 <- matrix(1) else Para$Sigma2 <- object$sigma2
  LambdaMean <- apply(object$lambda, 2, mean) #colMeans
  EtaMean <- apply(object$eta, 2, mean) #colMeans
  if (DatObj$P > 0) BetaMean <- apply(object$beta, 2, mean) else BetaMean <- matrix(0, ncol = 1, nrow = 0)
  if (!is.null(object$sigma2)) Sigma2Mean <- apply(object$sigma2, 2, mean) else Sigma2Mean <- matrix(1, nrow = 1, ncol = 1)
  Lambda <- matrix(LambdaMean, nrow = M * O, ncol = K, byrow = TRUE)
  Eta <- matrix(EtaMean, ncol = 1)
  MuMean <- array(kronecker(EyeNu, Lambda) %*% Eta + X %*% BetaMean, dim = c(M, O, Nu))
  Sigma2 <- t(matrix(Sigma2Mean, nrow = (O-C), ncol = M, byrow = TRUE))
  CovMean <- array(0, dim = c(M, O, Nu))
  count <- 1
  for (f in 1:O) {
    if (FamilyInd[f] %in% 0:2) {# if FamilyInd[f]==3, then the M*Nu entries CovMean[ , f, ] are null values (0 here)
      CovMean[ , f, ] <- matrix(Sigma2[, count], nrow = M, ncol = Nu)
      count <- count + 1
    }
  }
  Para$MuMean <- MuMean
  Para$CovMean <- CovMean
  
  ###Compute the Log-likelihood using our Rcpp function GetLogLik
  LogLik <- NULL
  if (("dic" %in% diags) | ("waic" %in% diags)) LogLik <- GetLogLik(DatObj, Para, NKeep, Verbose) #of length NKeep
  
  ###Compute DIC diagnostics
  dic <- NULL
  if ("dic" %in% diags) {
    LogLikMean <- GetLogLikMean(DatObj, Para) #scalar
    DBar <- -2 * mean(LogLik)
    DHat <- -2 * LogLikMean
    pD <- DBar - DHat
    DIC <- DBar + pD
    dic <- list(dic = DIC, pd = pD)    
  }
  
  ###Compute PPD diagnostics
  PPD <- dinfdiaglist <- meanIClist <- NULL
  if (("dinf" %in% diags) | ("meanIC" %in% diags)) {
    PPD <- SamplePPD(DatObj, Para, NKeep, Verbose) # of dim N x NKeep
    if ("dinf" %in% diags){
      PPDMean <- apply(PPD, 1, mean) # of length N = m x O x T (rowMeans)
      PPDVar <- apply(PPD, 1, var) # of length N 
      P <- sum(PPDVar)
      G <- sum( (PPDMean - YObserved) ^ 2) #lead to postMeanMSE
      DInf <- G + P
      dinfdiaglist <- list(p = P, g = G, dinf = DInf)
    }
    if ("meanIC" %in% diags){
      diffMat <- sweep(PPD, 1, YObserved, "-")
      meanIClist <- list(postMSE = mean(rowMeans(diffMat^2)), postVar = mean(apply(PPD, 1, var))) #postMSE and postVar
    }
  }
  
  ###Compute WAIC diagnostics
  waic <- NULL
  if ("waic" %in% diags) {
    lppd <- log(mean(exp(LogLik)))   
    if (is.infinite(lppd)) {
      M <- max(LogLik)
      lppd <- -log(dim(LogLik)[1]) + (M + log(sum(exp(LogLik - M))))
    }
    p_waic_1 <- 2 * (lppd - mean(LogLik) )
    p_waic_2 <- var(LogLik)
    waic <- -2 * lppd + 2 * p_waic_2
    waic <- list(waic = waic, p_waic = p_waic_2, lppd = lppd, p_waic_1 = p_waic_1)
  }
  
  if (!keepDeviance & !keepPPD) diags <- list(dic = dic, dinf = dinfdiaglist, meanIC = meanIClist, waic = waic)
  if (!keepDeviance & keepPPD) diags <- list(dic = dic, dinf = dinfdiaglist, meanIC = meanIClist, waic = waic, PPD = t(PPD))
  if (keepDeviance & !keepPPD) diags <- list(dic = dic, dinf = dinfdiaglist, meanIC = meanIClist, waic = waic, deviance = -2 * LogLik)
  if (keepDeviance & keepPPD) diags <- list(dic = dic, dinf = dinfdiaglist, meanIC = meanIClist, waic = waic, deviance = -2 * LogLik, PPD = t(PPD))
  return(diags)
}
