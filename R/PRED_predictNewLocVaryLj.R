###Function to get model fit diagnostics given a varyLjBFA object
#'
#' predictNewLocVaryLj
#'
#' Predicts observations at new location point(s) from a \code{\link{varyLjBFA}} or \code{\link{VAR1varyLjBFA}} model.
#'
#' @param object A \code{\link{varyLjBFA}} model object for which predictions
#'  are desired from.
#'
#' @param NNewLoc A positive integer indicating the number of new spatial locations for prediction.
#'
#' @param distOrigNew An \code{M x NNewLoc} dimensional distance matrix for distances between the original and new spatial points. 
#' If \code{include.space = FALSE}, then we don't need to specify \code{distOrigNew} and can leave it as its default value \code{NULL}. 
#' If \code{include.space = TRUE}, then this matrix must be specified no matter the value of the input \code{spatApprox}.
#' Since \code{spatial.structure} must be \code{"continuous"} to enable predictions at new spatial locations, the matrix should be a continuous distance matrix (e.g., Euclidean). 
#' 
#' @param distNewNew An \code{NNewLoc x NNewLoc} dimensional distance matrix for the new spatial points. 
#' If \code{include.space = FALSE} or \code{spatApprox = TRUE}, then we don't need to specify \code{distNewNew} and can leave it as its default value \code{NULL}. 
#' When \code{include.space = TRUE}, the matrix should be a continuous distance matrix (e.g., Euclidean), 
#' since \code{spatial.structure} must be \code{"continuous"} to enable predictions at new spatial locations.
#' 
#' @param NewX A matrix including covariates at times \code{1:T} for the new location points. 
#'  \code{NewX} must have dimension \code{(NNewLoc x O x T) x P}, where \code{NNewLoc} is the number of new location
#'  points being predicted. The default sets \code{NewX} to \code{NULL} and assumes that the covariates for all new locations
#'  are the same as the ones corresponding to the last reference location point.
#'  
#' @param NewTrials An array indicating the trials for categorical predictions. The array must have dimension \code{T x C x NNewLoc}
#'  and contain only non-negative integers. The default sets \code{NewTrials} to \code{NULL} and assumes that the trials for all predictions
#'  are the same as the ones corresponding to the final reference time point.
#'
#' @param Verbose A boolean logical indicating whether progress should be output.
#'
#' @param seed An integer value used to set the seed for the random number generator.
#'
#' @details \code{predictNewLocVaryLj} predicts vectors at new spatial location(s).
#'  The function returns the predicted factor loadings matrices \code{Lambda} and outcomes (\code{Y}).
#'  The function also returns the posterior predicted latent spatial vectors \code{alpha} and the corresponding \code{weights} for the
#'  new spatial location(s), which (when \code{include.space = TRUE}) are obtained via Bayesian kriging when \code{spatApprox = FALSE} 
#'  and via nearest-neighbor kriging when \code{spatApprox = TRUE}. 
#'
#' @return \code{predictNewLocVaryLj} returns a list containing the following objects.
#'
#'   \describe{
#'  
#'  \item{\code{Alpha}} {A \code{list} of posterior predicted values for \code{alpha}. The list contains \code{K x NKeep} elements (ordered by column), 
#'   each of which is a \code{matrix} of dimension \code{(NNewLoc x O) x (Lj - 1)} for \code{Lj} at that MCMC iteration. If \code{Lj = 1} for a certain \code{j}
#'   at a particular MCMC iteration, then that corresponding \code{matrix} in \code{Alpha} is of dimension \code{(NNewLoc x O) x 1} with entries all set to +Inf.
#'   In each element matrix for each clustering group, the corresponding entries are ordered first by observation type and then spatially.
#'   (the first O correspond to the first new location point for prediction, the next O correspond to the second new location point for prediction and so on)}
#'
#'   \item{\code{Weights}} {A \code{list} of posterior predicted values for \code{weights}. The list contains \code{K x NKeep} elements (ordered by column), 
#'   each of which is a \code{matrix} of dimension \code{(NNewLoc x O) x Lj} for \code{Lj} at that MCMC iteration.
#'   In each element matrix for each clustering group, the corresponding entries are ordered first by observation type and then spatially.
#'   (the first O correspond to the first new location point for prediction, the next O correspond to the second new location point for prediction and so on)}
#'
#'   \item{\code{Lambda}}{A \code{matrix} of dimension \code{NKeep x (NNewLoc x O x K)} containing the posterior predicted factor loadings matrices,
#'   where \code{O} is the number of observation types and \code{K} is the number of latent factors. 
#'   For each kept MCMC iteration, the corresponding predicted entries for the factor loadings matrix are ordered first by observation type,
#'   then spatially, and finally by factor. 
#'   (the first (NNewLoc x O) entries correspond to factor 1, the next (NNewLoc x O) entries correspond to factor 2 and so on;
#'   the first O rows correspond to the first new location point for prediction, the next O rows correspond to the second new location point for prediction and so on)}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{NNewLoc} posterior predictive distribution matrices.
#'   Each matrix is of dimension \code{NKeep x (T x O)}, where \code{T} is the number of time points and \code{O} is the number of observation types.
#'   For each kept MCMC iteration, the values are ordered first temporally and then by observation type (the first T correspond to the first observation type, 
#'   the next T correspond to the second observation type and so on). Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#' 
#' @references Yifan Cheng
#' @export

predictNewLocVaryLj <- function(object, NNewLoc, distOrigNew = NULL, distNewNew = NULL, NewX = NULL, NewTrials = NULL, Verbose = TRUE, seed = 27) {
  
  if (missing(object)) stop('"object" is missing')
  if (! (is.varyLjBFA(object) | is.VAR1varyLjBFA(object)) ) stop('"object" must be of class varyLjBFA or VAR1varyLjBFA')
  if (!is.logical(Verbose)) stop('"Verbose" must be a logical')
  if (missing(NNewLoc)) stop("NNewLoc: missing")
  if (!is.scalar(NNewLoc)) stop('NNewLoc must be a scalar')
  if (is.na(NNewLoc)) stop('NNewLoc should not be NA')
  if (is.infinite(NNewLoc)) stop('NNewLoc should not be infinite')
  if (!is.wholenumber(NNewLoc) | NNewLoc <= 0) stop('NNewLoc must be a strictly positive integer')

  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  K <- DatObj$K
  C <- DatObj$C
  h <- DatObj$h
  FamilyInd <- DatObj$FamilyInd 
  IS <- DatObj$IS
  SA <- DatObj$SA
  SpCorInd <- DatObj$SpCorInd
  
  if (DatObj$spatPred == 0) stop("'storeSpatPredPara' must be TRUE to enable predictions at new spatial locations")
  if (IS == 1) {
    if (SpCorInd == 1) stop("If a spatial process is included when fitting the model, then 'spatial.structure' must be 'continuous' to enable predictions at new spatial locations")
    if (is.null(distOrigNew)) stop("When include.space = TRUE, 'distOrigNew' must be specified no matter the value of the input 'spatApprox'")
    if ((SA == 0) & is.null(distNewNew)) stop("When include.space = TRUE and spatApprox = FALSE, 'distNewNew' must be specified")
  } 
  if (!is.null(distOrigNew)){
    if (!is.matrix(distOrigNew)) stop("'distOrigNew' must be a matrix")
    if (nrow(distOrigNew) != M) stop(paste0("'distOrigNew' must be an M x NNewLoc dimensional matrix"))
    if (ncol(distOrigNew) != NNewLoc) stop("'distOrigNew' must be an M x NNewLoc dimensional matrix")
    if (any(is.na(distOrigNew))) stop("'distOrigNew' should not contain missing values")
    if (any(is.infinite(distOrigNew))) stop("'distOrigNew' should not contain infinite values")
    if (any(distOrigNew < 0)) stop("'distOrigNew' must not contain negative entries")
  }
  if (!is.null(distNewNew)){
    if (!is.matrix(distNewNew)) stop("'distNewNew' must be a matrix")
    if (nrow(distNewNew) != NNewLoc) stop(paste0("'distNewNew' must be an NNewLoc x NNewLoc dimensional matrix"))
    if (ncol(distNewNew) != NNewLoc) stop("'distNewNew' must be an NNewLoc x NNewLoc dimensional matrix")
    if (any(is.na(distNewNew))) stop("'distNewNew' should not contain missing values")
    if (any(is.infinite(distNewNew))) stop("'distNewNew' should not contain infinite values")
    if (any(distNewNew < 0)) stop("'distNewNew' must not contain negative entries")
    if (any(distNewNew != t(distNewNew))) stop('distNewNew must be symmetric')
  }
  
  if (P > 0){
    spatIndicesSub <- rep(1:M, O)
    spatIndices <- rep(spatIndicesSub, Nu)
    if (is.null(NewX)) {
      XM <- object$datobj$X[spatIndices == M, ]
      XMresize <- matrix(XM, O, Nu * P)
      NewXresize <- do.call("rbind", rep(list(XMresize), NNewLoc))#the first T x O obs from new location point 1; the first O from time point 1
      NewX <- matrix(NewXresize, (NNewLoc * O * Nu), P) #the first (NNewLoc x O) obs from time point 1; the first O from new loc 1
    } else {
      if (!is.matrix(NewX)) stop('NewX must be a matrix')
      if (nrow(NewX) != (NNewLoc * O * Nu)) stop("NewX: Must be a matrix with dimension (NNewLoc x O x T) x P")
      if (ncol(NewX) != P) stop("NewX: Must be a matrix with dimension (NNewLoc x O x T) x P")
      if (any(is.na(NewX))) stop("NewX should not contain missing values")
      if (any(is.infinite(NewX))) stop("NewX must contain strictly finite entries")
      NewX <- NewX
    }
  } else { # when P = 0
    NewX <- matrix(0, ncol = 0, nrow = (NNewLoc * O * Nu))
  }
  
  if (C == 0) {
    DatObj$Trials <- array(0, dim = c(Nu, O, NNewLoc))
  } else if ((C > 0) & is.null(NewTrials)) {
    Trials <- array(dim = c(Nu, C, NNewLoc))
    for (n in 1:NNewLoc) Trials[ , , n] <- t(DatObj$Trials[M, which(FamilyInd == 3), ]) 
    # the original datobj DatObj.Trials is of dimension M*O*Nu if C>0 with entries corresponding to non-count observation types set to 1 (see MCMC_Create.R) 
    DatObj$Trials <- Trials
  } else {#if (C > 0) & !is.null(NewTrials) 
    Trials <- NewTrials
    if (!is.array(Trials)) stop('The NewTrials input must be an array if it is not NULL')
    if (dim(Trials)[1] != Nu) stop("NewTrials: Must be an array with dimension T x C x NNewLoc")
    if (dim(Trials)[2] != C) stop("NewTrials: Must be an array with dimension T x C x NNewLoc")
    if (dim(Trials)[3] != NNewLoc) stop("NewTrials: Must be an array with dimension T x C x NNewLoc")
    if (any(is.na(Trials))) stop("Trials may contain no missing values")
    if (any(!is.finite(Trials))) stop("Trials must contain strictly finite entries")
    if (any(Trials != floor(Trials))) stop("Trials must contain integers only")
    if (any(Trials < 0)) stop("Trials must contain non-negative integers only")
    DatObj$Trials <- Trials
  }
  
  nnIndpred <- matrix(0, NNewLoc, h)
  if ((IS == 1) & (SA == 1)){
    for(i in 1:NNewLoc){
      nnIndpred[i, ] <- sort(order(distOrigNew[, i])[1:h])
    }
  }
  
  if(is.null(distOrigNew)) distOrigNew <- matrix(1,1,1)
  if(is.null(distNewNew)) distNewNew <- matrix(1,1,1)
  NKeep <- dim(object$lambda)[1]
  HyPara <- object$hypara
  DatObj$sigma2hyparaA <- HyPara$A
  DatObj$sigma2hyparaB <- HyPara$B
  DatObj$NewX <- NewX
  DatObj$NNewLoc <- NNewLoc
  DatObj$distOrigNew <- distOrigNew
  DatObj$distNewNew <- distNewNew
  ##DatObj$storeW <- 1 # weights are always stored for the bfaVaryingLjs main function
  DatObj$nnIndpred = nnIndpred - 1
  
  Para <- list()
  Para$Eta <- object$eta
  if (DatObj$P > 0) Para$Beta <- object$beta else Para$Beta <- matrix(0, ncol = 0, nrow = NKeep)
  Para$Rho <- object$rho
  Para$Kappa <- object$kappa
  Para$Lambda <- object$lambda
  Para$Theta <- object$theta
  Para$Alpha <- object$alpha
  Para$LjVec <- object$ljvec

  set.seed(seed)
  LatentKrig <- AlphaKriggingVaryLj(DatObj, Para, NKeep, Verbose)
  AlphaKrig <- LatentKrig$alpha
  if (DatObj$storeW == 1) {WeightsKrig <- LatentKrig$weights} else {WeightsKrig <- NULL}
  LambdaKrig <- LatentKrig$lambda

  YKrig <- YKriggingSpat(DatObj, Para, LambdaKrig, NKeep, Verbose)
  YOut <- list()
  YInd <- expand.grid(1:Nu, 1:O)
  for (n in 1:NNewLoc) {
    YOut[[n]] <- t(YKrig[, n, ])
    colnames(YOut[[n]]) <- paste("Y", M + n, YInd[, 2], YInd[, 1], sep="_") # the first T corresponds to obs type 1 and so on
  }
  names(YOut) <- paste0("Y", M + 1:NNewLoc)  
  return(list(Alpha = AlphaKrig, Weights = WeightsKrig, Lambda = t(LambdaKrig), Y = YOut))
  
}





### Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
