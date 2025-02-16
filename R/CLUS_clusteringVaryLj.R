###Function to obtain spatial temporal trends clustering results given a varyLjBFA object
#'
#' clusteringVaryLj
#'
#' Cluster spatial points into regions with similar temporal trajectories from a \code{\link{varyLjBFA}} model.
#'
#' @param object A \code{\link{varyLjBFA}} or \code{\link{VAR1varyLjBFA}} model object on which sptial temporal trends clustering is based.
#' 
#' @param o A positive integer less than or equal to O (the number of observation types) indicating the single observation type that we want to perform clustering for.
#'
#' @param nkeep A positive integer less than or equal to NKeep (the number of post-burn-in post-thinned MCMC iterations) indicating the number of selected kept 
#' post-burn-in MCMC iterations whose posterior parameter estimates our spatial temporal trends clustering will be based on. 
#' 
#' @param nCent A positive integer indicating the number of clusters for k-means clustering.
#'
#' @return \code{clusteringVaryLj} returns an object from R's \code{kmeans()} function, which is a list with elements \code{cluster}, \code{centers},
#' \code{totss}, \code{withinss}, \code{tot.withinss}, \code{betweenss}, \code{size}, \code{iter}, and \code{ifault}.
#' 
#' @references Yifan Cheng
#' @export

clusteringVaryLj <- function(object, o = 1, nkeep = 1, nCent = 3){
  
  if (missing(object)) stop('"object" is missing')
  if (! (is.varyLjBFA(object) | is.VAR1varyLjBFA(object)) ) stop('"object" must be of class varyLjBFA or VAR1varyLjBFA')
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
  if (!is.wholenumber(nCent) | nCent < 0) stop('"nCent" must be a positive integer')
  
  NKeep <- dim(object$lambda)[1]
  DatObj <- object$datobj
  M <- DatObj$M
  O <- DatObj$O
  K <- DatObj$K
  selectedIter <- seq(from = 1, by = floor((NKeep-1)/nkeep)+1, length.out = nkeep)
  
  if (o > O) stop('o must be no larger than O')
  if (nkeep > NKeep) stop('nkeep must be no larger than NKeep')
  if (DatObj$storeW == 0) stop('"storeWeights" must be TRUE if we want to perform temporal trends clustering from the output object')
  
  postWeightsAll <- object$weights # a list that contains \code{K x NKeep} elements (ordered by column), 
  # each of which is a matrix of dimension Lj x (M x O) for Lj at that MCMC iteration.
  uppervec <- selectedIter * K
  lowervec <- (selectedIter - 1) * K + 1
  keptInd <- as.vector(mapply(`:`, lowervec, uppervec))
  postWeightsList <- postWeightsAll[keptInd] # a sublist containing only element matrices corresponding to the selected MCMC iterations
  postWeightsVec <- unlist(lapply(postWeightsList, t))
  postWeightsArray <- array(postWeightsVec, dim = c(M, O, length(postWeightsVec)/(M*O)))
  postWeightsMat <- matrix(postWeightsArray[, o, ], nrow = M) # of dim M x (sum of the Lj's for all j and all selected iter)
  
  return(kmeans(postWeightsMat, centers = nCent)) # returns a kmeans(w) object
  
}



### Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
