###Function to get model fit diagnostics given a FixedLbfaVAR1 or VAR1varyLjBFA object
#'
#' predictNewTime
#'
#' Predicts future observation(s) from a \code{\link{FixedLbfaVAR1}} or \code{\link{VAR1varyLjBFA}} model.
#'
#' @param object A \code{\link{FixedLbfaVAR1}} or \code{\link{VAR1varyLjBFA}} model object for which predictions
#'  are desired from.
#'
#' @param NNewTime An integer indicating the desired number of consecutive future time point(s) after T for prediction. 
#' 
#' @param NewX A matrix including covariates at times \code{NewTimes} for prediction. 
#'  \code{NewX} must have dimension \code{(M x O x NNewTime) x P}. Where \code{NNewTime} is the number of temporal 
#'  locations being predicted. The default sets \code{NewX} to \code{NULL}, which assumes that the covariates for all predictions 
#'  are the same as the final time point.
#'  
#' @param NewTrials An array indicating the trials for categorical predictions. The array must have dimension \code{M x C x NNewTime}
#'  and contain only non-negative integers. The default sets \code{NewTrials} to \code{NULL}, which assumes that trials for all predictions
#'  are the same as the final time point.
#'
#' @param Verbose A boolean logical indicating whether progress should be output.
#'
#' @param seed An integer value used to set the seed for the random number generator
#'  (default = 54).
#'  
#' @param ... other arguments.
#'
#' @details \code{predictNewTime} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged factors (\code{Eta}) and also the observed outcomes (\code{Y}).
#'
#' @return \code{predictNewTime} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{Eta}}{A \code{list} containing \code{NNewTime} matrices, one for each new time prediction. Each matrix is dimension \code{NKeep x K}, where
#'   \code{K} is the number of latent factors. Each matrix contains posterior samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{NNewTime} posterior predictive distribution matrices. 
#'   Each matrix is of dimension \code{NKeep x (M * O)}, where \code{M} is the number of spatial locations and \code{O} the number of observation types.
#'   Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#' 
#' @references Yifan Cheng
#' @export

VAR1predictNewTime <- function(object, NNewTime, NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 54, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!(is.FixedLbfaVAR1(object) | is.VAR1varyLjBFA(object)) ) stop('"object" must be of class FixedLbfaVAR1 or VAR1varyLjBFA')
  if (missing(NNewTime)) stop('"NNewTime" is missing')
  if (!is.scalar(NNewTime)) stop('"NNewTime" must be a scalar')
  if (is.na(NNewTime)) stop('"NNewTime" cannot be NA')
  if (!is.finite(NNewTime)) stop('"NNewTime" cannot be infinite')
  if (!is.wholenumber(NNewTime) | NNewTime <= 0) stop('"NNewTime" must be a strictly positive integer')
  if (!is.logical(Verbose)) stop('"Verbose" must be a logical')
  
  ###Set seed for reproducibility
  set.seed(seed)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  K <- DatObj$K
  FamilyInd <- DatObj$FamilyInd 

  ###Get covariates
  if (is.null(NewX)) {
    XNu <- object$datobj$X[object$datobj$Indeces == (object$datobj$Nu - 1), ]
    NewX <- do.call("rbind", rep(list(XNu), NNewTime))
  } else {
    if (!is.matrix(NewX)) stop('NewX must be a matrix')
    if (nrow(NewX) != (M * O * NNewTime)) stop("NewX: Must be a matrix with dimension (M x O x NNewTime) x P")
    if (ncol(NewX) != P) stop("NewX: Must be a matrix with dimension (M x O x NNewTime) x P")
    if (any(is.na(NewX))) stop("NewX may contain no missing values")
    if (any(is.infinite(NewX))) stop("NewX must contain strictly finite entries")
    NewX <- NewX
  }
  
  ###Update DatObj
  DatObj$NNewTime <- NNewTime
  DatObj$NewX <- NewX
  
  ###Create Trials object
  if (DatObj$C == 0) {
    DatObj$Trials <- array(0, dim = c(M, O, NNewTime))
  } else if ((DatObj$C > 0) & is.null(NewTrials)) {
    Trials <- array(dim = c(M, DatObj$C, NNewTime))
    for (n in 1:NNewTime) Trials[, , n] <- DatObj$Trials[, which(FamilyInd == 3), Nu] 
    # the original datobj DatObj. Trials is of dimension M*O*Nu if C>0 with entries corresponding to non-count observation types set to 1 (see MCMC_Create.R) 
    DatObj$Trials <- Trials
  } else {# if (DatObj$C > 0) & !is.null(NewTrials)
    Trials <- NewTrials
    if (!is.array(Trials)) stop('The NewTrials input must be an array if it is not NULL')
    if (dim(Trials)[1] != M) stop("NewTrials: Must be an array with dimension M x C x NNewTime")
    if (dim(Trials)[2] != DatObj$C) stop("NewTrials: Must be an array with dimension M x C x NNewTime")
    if (dim(Trials)[3] != NNewTime) stop("NewTrials: Must be an array with dimension M x C x NNewTime")
    if (any(is.na(Trials))) stop("Trials may contain no missing values")
    if (any(is.infinite(Trials))) stop("Trials must contain strictly finite entries")
    if (any(Trials != floor(Trials))) stop("Trials must contain integers only")
    if (any(Trials < 0)) stop("Trials must contain non-negative integers only")
    DatObj$Trials <- Trials
  }
  
  ###Set mcmc object
  NKeep <- dim(object$lambda)[1]

  ###Create parameter object
  Para <- list()
  Para$A <- object$A
  Para$Upsilon <- object$upsilon
  Para$Lambda <- object$lambda
  Para$Eta <- object$eta
  if (DatObj$P > 0) Para$Beta <- object$beta else Para$Beta <- matrix(0, ncol = 0, nrow = NKeep)
  if (is.null(object$sigma2)) Para$Sigma2 <- matrix(1) else Para$Sigma2 <- object$sigma2
  
  ###Obtain samples of eta using Bayesian krigging
  EtaKrig <- VAR1EtaKrigging(DatObj, Para, NKeep, Verbose)

  ###Obtain samples of observed Y
  YKrig <- VAR1YKriggingTemp(DatObj, Para, EtaKrig, NKeep, Verbose)

  ###Format eta samples for output
  EtaOut <- list()
  Eta <- array(t(EtaKrig), dim = c(NKeep, K, NNewTime))
  for (n in 1:NNewTime) {
    EtaOut[[n]] <- Eta[, , n]
    colnames(EtaOut[[n]]) <- paste0("Eta", 1:K, "_", Nu + n)
  }
  names(EtaOut) <- paste0("Eta", Nu + 1:NNewTime)
  
  ###Format Y samples for output
  YOut <- list()
  YInd <- expand.grid(1:M, 1:O)
  for (n in 1:NNewTime) {
    YOut[[n]] <- t(YKrig[, n, ])
    colnames(YOut[[n]]) <- paste("Y", Nu + n, YInd[, 2], YInd[, 1], sep="_")
  }
  names(YOut) <- paste0("Y", Nu + 1:NNewTime)  
    
  ###Return formated samples
  return(list(Eta = EtaOut, Y = YOut))

}
