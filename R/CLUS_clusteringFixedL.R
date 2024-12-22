###Function to obtain spatial temporal trends clustering results given a FixedLbfa object
#'
#' clusteringFixedL
#'
#' Cluster spatial points into regions with similar temporal trajectories from a \code{\link{FixedLbfa}} model.
#'
#' @param object A \code{\link{FixedLbfa}} or \code{\link{FixedLbfaVAR1}} model object on which sptial temporal trends clustering is based.
#' 
#' @param o A positive integer less than or equal to O (the number of observation types) indicating the single observation type that we want to perform clustering for.
#'
#' @param nkeep A positive integer less than or equal to NKeep (the number of post-burn-in post-thinned MCMC iterations) indicating the number of selected kept 
#' post-burn-in MCMC iterations whose posterior parameter estimates our spatial temporal trends clustering will be based on. If \code{nkeep = 0}, the posterior
#' mean weights estimates will be used for our clustering analysis.
#' 
#' @param nCent A positive integer indicating the number of clusters for k-means clustering.
#'
#' @return \code{clusteringFixedL} returns an object from R's \code{kmeans()} function, which is a list with elements \code{cluster}, \code{centers},
#' \code{totss}, \code{withinss}, \code{tot.withinss}, \code{betweenss}, \code{size}, \code{iter}, and \code{ifault}.
#' 
#' @references Yifan Cheng
#' @export


clusteringFixedL <- function(object, o = 1, nkeep = 1, nCent = 3) {
  
  if (missing(object)) stop('"object" is missing')
  if (! (is.FixedLbfa(object) | is.FixedLbfaVAR1(object)) ) stop('"object" must be of class FixedLbfa or FixedLbfaVAR1')
  if (!is.scalar(o)) stop('o must be a scalar')
  if (is.na(o)) stop('o cannot be NA')
  if (!is.finite(o)) stop('o cannot be infinite')
  if (!is.wholenumber(o) | o < 0) stop('o must be a positive integer')
  if (!is.scalar(nkeep)) stop('nkeep must be a scalar')
  if (is.na(nkeep)) stop('nkeep cannot be NA')
  if (!is.finite(nkeep)) stop('nkeep cannot be infinite')
  if (!is.wholenumber(nkeep) | nkeep < 0) stop('nkeep must be a positive integer')
  if (!is.scalar(nCent)) stop('nCent must be a scalar')
  if (is.na(nCent)) stop('nCent cannot be NA')
  if (!is.finite(nCent)) stop('nCent cannot be infinite')
  if (!is.wholenumber(nCent) | nCent < 0) stop('nCent must be a positive integer')
  
  NKeep <- dim(object$lambda)[1]
  DatObj <- object$datobj
  M <- DatObj$M
  O <- DatObj$O
  K <- DatObj$K
  L <- DatObj$L
  selectedIter <- seq(from = 1, by = floor((NKeep-1)/nkeep)+1, length.out = nkeep)

  if (o > O) stop('o must be no larger than O')
  if (nkeep > NKeep) stop('nkeep must be no larger than NKeep')
  if (DatObj$CL == 0) stop('"clustering" must be TRUE if we want to perform temporal trends clustering from the output object')
  if (DatObj$storeW == 0) stop('"storeWeights" must be TRUE if we want to perform temporal trends clustering from the output object')
  
  if (DatObj$alphaWsFiles == 0) {
    postWeightsAll <- object$weights # of dim NKeep x (M x O x K x L)
    if (nkeep > 0){
      postWeightsSelectedIter <- postWeightsAll[selectedIter,] # of dim nkeep x (M x O x K x L)
      postWeightsSelectedIterArrayO <- array(postWeightsSelectedIter, dim = c(nkeep*M, O, K*L))[, o, ] 
      postWeightsArray <- aperm(array(postWeightsSelectedIterArrayO, dim = c(nkeep, M, K*L)), c(2, 1, 3)) # M x nkeep x (K*L)
      postWeightsMat <- matrix(postWeightsArray, M, nkeep*K*L)
    } else{ # if nkeep == 0
      postMeanWeights <- colMeans(postWeightsAll) # of dim 1 x (M x O x K x L)
      postMeanWeightsArrayO <- array(postMeanWeights, dim = c(M, O, K*L))[, o, ] 
      postWeightsMat <- matrix(postMeanWeightsArrayO, M, K*L)
    }
  } else {
    if (nkeep > 0){
      postWeightsList <- list()
      for(n in 1:nkeep){
        matrixn <- t(as.matrix(read.table(paste("fixedLweightsKeptIter", selectedIter[n], ".txt", sep=""), skip = 2))) # of dim (M*O) x (K*L)
        postWeightsList[[n]] <- matrix(array(matrixn, dim = c(M, O, K*L))[, o, ], nrow = M) # of dim M x (KL)
      } #read in all posterior weights text files for selected kept iterations
      postWeightsMat <- matrix(unlist(postWeightsList), nrow = M) # of dim M x (nkeep x K x L)
    } else{  # if nkeep == 0
      postWeightsList <- list()
      for(n in 1:NKeep){
        matrixn <- t(as.matrix(read.table(paste("fixedLweightsKeptIter", selectedIter[n], ".txt", sep=""), skip = 2))) # of dim (M*O) x (K*L)
        postWeightsList[[n]] <- matrix(array(matrixn, dim = c(M, O, K*L))[, o, ], nrow = M) # of dim M x (KL)
      } #read in all posterior weights text files 
      postWeightsAllIter <- matrix(unlist(postWeightsList), nrow = M * K * L) # of dim (M x K x L) x NKeep
      postWeightsMat <- matrix(rowMeans(postWeightsAllIter), nrow = M) # of dim M x (K x L)
    }
  } 
  
  return(kmeans(postWeightsMat, centers = nCent)) # returns a kmeans(w) object
  
}



### Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
