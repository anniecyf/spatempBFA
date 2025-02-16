###Function to get model fit diagnostics given a FixedLbfa or varyLjBFA object
#'
#' predictNewTime
#'
#' Predicts future observation(s) from a \code{\link{FixedLbfa}} or \code{\link{varyLjBFA}} model.
#'
#' @param object A \code{\link{FixedLbfa}} or \code{\link{varyLjBFA}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes A numeric vector including desired time point(s) for prediction. If \code{include.time = equalTimeDist = TRUE} when getting this function's input \code{object} from the main function \code{bfa_sp} 
#' and the pre-normalized time distance does not equal 1, then we should standardize the new time points accordingly. \code{NewTimes} should be the standardized desired time point(s) for prediction.
#' 
#' @param NewX A matrix including covariates at times \code{NewTimes} for prediction. 
#'  \code{NewX} must have dimension \code{(M x O x NNewVisits) x P}. Where \code{NNewVisits} is the number of temporal 
#'  locations being predicted. The default sets \code{NewX} to \code{NULL}, which assumes that the covariates for all predictions 
#'  are the same as the final time point.
#'  
#' @param NewTrials An array indicating the trials for categorical predictions. The array must have dimension \code{M x C x NNewVisits}
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
#'   \item{\code{Eta}}{A \code{list} containing \code{NNewVisits} matrices, one for each new time prediction. Each matrix is dimension \code{NKeep x K}, where
#'   \code{K} is the number of latent factors. Each matrix contains posterior samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{NNewVisits} posterior predictive distribution matrices. 
#'   Each matrix is of dimension \code{NKeep x (M * O)}, where \code{M} is the number of spatial locations and \code{O} the number of observation types.
#'   Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#' 
#' @references Yifan Cheng
#' @references Berchuck, S. I., Janko, M., Medeiros, F. A., Pan, W., & Mukherjee, S. (2021). Bayesian Non-Parametric Factor Analysis for Longitudinal Spatial Surfaces. Bayesian Analysis, 17(2), 1–30.
#' @export

predictNewTime <- function(object, NewTimes, NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 54, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!(is.FixedLbfa(object) | is.varyLjBFA(object))) stop('"object" must be of class FixedLbfa or varyLjBFA')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may contain no missing values")
  if (any(is.infinite(NewTimes))) stop("NewTimes must contain strictly finite entries")
  if (any(NewTimes < 0)) stop('NewTimes vector contains at least one negative entry')
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

  ###Create updated distance matrix
  TimeFixed <- DatObj$Time
  Time <- sort(c(TimeFixed, NewTimes))
  TimeDist <- abs(outer(Time, Time, "-" ))
  NNewVisits <- length(NewTimes)
  NewVisits <- OriginalVisits <- NULL
  for (i in 1:NNewVisits) NewVisits <- c(NewVisits, which(NewTimes[i] == Time) - 1)
  for (i in 1:Nu) OriginalVisits <- c(OriginalVisits, which(TimeFixed[i] == Time) - 1)

  ###Get covariates
  if (is.null(NewX)) {
    XNu <- object$datobj$X[object$datobj$Indeces == (object$datobj$Nu - 1), ]
    NewX <- do.call("rbind", rep(list(XNu), NNewVisits))
  } else {
    if (!is.matrix(NewX)) stop('NewX must be a matrix')
    if (nrow(NewX) != (M * O * NNewVisits)) stop("NewX: Must be a matrix with dimension (M x O x NNewVisits) x P")
    if (ncol(NewX) != P) stop("NewX: Must be a matrix with dimension (M x O x NNewVisits) x P")
    if (any(is.na(NewX))) stop("NewX may contain no missing values")
    if (any(is.infinite(NewX))) stop("NewX must contain strictly finite entries")
    NewX <- NewX
  }
  
  ###Update DatObj
  DatObj$NewVisits <- NewVisits
  DatObj$OriginalVisits <- OriginalVisits
  DatObj$TimeDist <- TimeDist
  DatObj$NNewVisits <- NNewVisits
  DatObj$EyeK <- diag(DatObj$K)
  DatObj$NewX <- NewX
  
  ###Create Trials object
  if (DatObj$C == 0) {
    DatObj$Trials <- array(0, dim = c(M, O, NNewVisits))
  } else if ((DatObj$C > 0) & is.null(NewTrials)) {
    Trials <- array(dim = c(M, DatObj$C, NNewVisits))
    for (n in 1:NNewVisits) Trials[, , n] <- DatObj$Trials[, which(FamilyInd == 3), Nu] 
    # the original datobj DatObj. Trials is of dimension M*O*Nu if C>0 with entries corresponding to non-count observation types set to 1 (see MCMC_Create.R) 
    DatObj$Trials <- Trials
  } else {# if (DatObj$C > 0) & !is.null(NewTrials)
    Trials <- NewTrials
    if (!is.array(Trials)) stop('The NewTrials input must be an array if it is not NULL')
    if (dim(Trials)[1] != M) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (dim(Trials)[2] != DatObj$C) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (dim(Trials)[3] != NNewVisits) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
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
  Para$Psi <- object$psi
  Para$Upsilon <- object$upsilon
  Para$Lambda <- object$lambda
  Para$Eta <- object$eta
  if (DatObj$P > 0) Para$Beta <- object$beta else Para$Beta <- matrix(0, ncol = 0, nrow = NKeep)
  if (is.null(object$sigma2)) Para$Sigma2 <- matrix(1) else Para$Sigma2 <- object$sigma2
  
  ###Obtain samples of eta using Bayesian krigging
  EtaKrig <- EtaKrigging(DatObj, Para, NKeep, Verbose)

  ###Obtain samples of observed Y
  YKrig <- YKriggingTemp(DatObj, Para, EtaKrig, NKeep, Verbose)

  ###Format eta samples for output
  EtaOut <- list()
  Eta <- array(t(EtaKrig), dim = c(NKeep, K, NNewVisits))
  for (n in 1:NNewVisits) {
    EtaOut[[n]] <- Eta[, , n]
    colnames(EtaOut[[n]]) <- paste0("Eta", 1:K, "_", Nu + n)
  }
  names(EtaOut) <- paste0("Eta", Nu + 1:NNewVisits)
  
  ###Format Y samples for output
  YOut <- list()
  YInd <- expand.grid(1:M, 1:O)
  for (n in 1:NNewVisits) {
    YOut[[n]] <- t(YKrig[, n, ])
    colnames(YOut[[n]]) <- paste("Y", Nu + n, YInd[, 2], YInd[, 1], sep="_")
  }
  names(YOut) <- paste0("Y", Nu + 1:NNewVisits)  
    
  ###Return formated samples
  return(list(Eta = EtaOut, Y = YOut))

}
