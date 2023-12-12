#' Spatial factor analysis using a Bayesian hierarchical model.
#'
#' \code{VAR1bfaFixedL} is a Markov chain Monte Carlo (MCMC) sampler for a Bayesian spatial factor analysis model. The spatial component is 
#' introduced using a Probit stick-breaking process prior on the factor loadings. The model is implemented using a Bayesian hierarchical framework.
#'
#' @param formula A \code{formula} object, corresponding to the spatial factor analysis model. The response must be on the left of a \code{~} operator, and the terms on the right 
#'                 must indicate the covariates to be included in the fixed effects. If no covariates are desired a zero should be used, \code{~ 0}. 
#'                 
#' @param data A required \code{data.frame} containing the variables (\code{Y}, additional \code{x} coviariate(s) if there is/are, and a \code{trials} variable if there
#'  is/are observation type(s) from \code{family "binomial"}) in the model. The data frame must contain \code{M x O x Nu} rows.
#'  Here, \code{M} represents the number of spatial locations, \code{O} the number of different observation types
#'  and \code{Nu} the number of temporal visits. The observations must be first be
#'  ordered spatially, second by observation type and then temporally. This means that the first \code{M x O} observations come from the first time point and
#'  the first \code{M} observations come the first spatial observation type.
#'  If there is/are observation type(s) from \code{family "binomial"}, then \code{trials} contains the numbers of trials as positive integers for each of the binomial observations.
#'  Entries in \code{trials} corresponding to non-binomial data can be specified as any arbitrary positive integer. The function will change these values to 1.
#'
#' @param family A character string or a vector of length \code{O} (if \code{O > 1}) of character strings indicating the distribution(s) of the observed data. Options for each observation type
#'  include: \code{"normal"}, \code{"probit"}, \code{"tobit"}, and \code{"binomial"}. If \code{O > 1} and \code{family} is of length 1, then all of the \code{O} observation types are from the same 
#'  \code{family} (distribution). Any combination of likelihoods can be used.
#'  
#' @param dist A \code{M x M} dimensional distance matrix. For a \code{discrete} spatial process the matrix contains binary adjacencies that dictate the
#'  spatial neighborhood structure and for \code{continuous} spatial processes the matrix should be a continuous distance matrix (e.g., Euclidean).
#'  
#' @param Nu A positive integer representing the number of evenly dispersed time points corresponding to our observed data.
#'
#' @param K A scalar that indicates the dimension (i.e., quantity) of latent factors.
#'  
#' @param L A positive integer indicating the number of latent clusters for each column of the factor loadings matrix. 
#'  
#' @param starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon}, 
#'  or \code{A} containing appropriate objects. \code{Beta} (or \code{Delta}) must either be a \code{P} (or \code{K}) dimensional
#'  vector or a scalar (the scalar populates the entire vector). \code{Sigma2} must be either a \code{M x (O - C)} matrix or a scalar.
#'  \code{Kappa} must be a \code{O x O} dimensional matrix, \code{Rho} a scalar, and \code{Upsilon} a \code{K x K} matrix.
#'
#' @param hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in (Berchuck et al.).
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, or \code{Upsilon} 
#'  containing further hyperparameter information. These objects are themselves
#'  \code{lists} and may be constructed as follows.
#'
#'  \code{Beta} is a \code{list} with two objects, \code{MuBeta} and \code{SigmaBeta}. These values represent the prior mean and variance 
#'  parameters for the multivariate normal prior.
#'  
#'  \code{Delta} is a \code{list} with two objects, \code{A1} and \code{A2}. These values represent the prior shape 
#'  parameters for the multiplicative Gamma shrinkage prior.
#'  
#'  \code{Sigma2} is a \code{list} with two objects, \code{A} and \code{B}. These values represent the shape and scale for the variance parameters.
#'  
#'  \code{Kappa} is a \code{list} with two objects,
#'  \code{SmallUpsilon} and \code{BigTheta}. \code{SmallUpsilon} represents the degrees of freedom parameter for the
#'  inverse-Wishart hyperprior and must be a real number scalar, while \code{BigTheta} represents
#'  the scale matrix and must be a \code{O x O} dimensional positive definite matrix.
#'  
#'  \code{Rho} is a \code{list} with two objects, \code{ARho} and \code{BRho}. \code{ARho}
#'  represents the lower bound for the uniform hyperprior, while \code{BRho} represents
#'  the upper bound. The bounds must be specified carefully. This is only specified for continuous spatial processes.
#'  
#'  \code{Upsilon} is a \code{list} with two objects,
#'  \code{Zeta} and \code{Omega}. \code{Zeta} represents the degrees of freedom parameter for the
#'  inverse-Wishart hyperprior and must be a real number scalar, while \code{Omega} represents
#'  the scale matrix and must be a \code{K x K} dimensional positive definite matrix.
#'  
#' @param tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the MCMC Metropolis steps. If \code{NULL} is not chosen then all
#'  of the tuning values must be specified.
#'
#'  When \code{NULL} is chosen then default tuning values are automatically generated to
#'  \code{1}. Otherwise a \code{list} must be provided with name \code{Rho}. 
#'  Each of these entries must be scalars containing tuning variances for their corresponding Metropolis updates.
#'
#' @param mcmc Either \code{NULL} or a \code{list} containing input values to be used
#'  for implementing the MCMC sampler. If \code{NULL} is not chosen then all
#'  of the MCMC input values must be specified.
#'
#'  \code{NBurn}: The number of sampler scans included in the burn-in phase. (default =
#'  \code{10,000})
#'
#'  \code{NSims}: The number of post-burn-in scans for which to perform the
#'   sampler. (default = \code{10,000})
#'
#'  \code{NThin}: Value such that during the post-burn-in phase, only every
#'  \code{NThin}-th scan is recorded for use in posterior inference (For return values
#'  we define, NKeep = NSims / NThin (default = \code{1}).
#'
#'  \code{NPilot}: The number of times during the burn-in phase that pilot adaptation
#'  is performed (default = \code{20})
#'
#' @param spatial.structure Character string indicating the type of spatial process. Options include:
#'  \code{"continuous"} (i.e., Gaussian process with exponential kernel) and \code{"discrete"} (i.e., proper CAR).
#'
#' @param seed An integer value used to set the seed for the random number generator
#'  
#' @param gamma.shrinkage A logical indicating whether a gamma shrinkage process prior is used for the variances of the factor loadings columns. If FALSE,
#'  the hyperparameters (A1 and A2) indicate the shape and rate for a gamma prior on the precisions. Default is TRUE. It can only be TRUE when \code{clustering = TRUE}.
#'
#' @param include.space A logical indicating whether a spatial process should be included. Default is TRUE, however if FALSE the spatial correlation matrix 
#'  \code{F(\rho)} is fixed as the \code{M x M} identity matrix. This specification overrides the inputs \code{spatial.structure}, \code{spatApprox}, and \code{alphaMethod}.
#'  
#' @param clustering A logical indicating whether the Bayesian non-parametric process should be used, default is TRUE. If FALSE is specified
#'  each column is instead modeled with an independent spatial process.
#'  
#' @param spatApprox A logical indicating whether spatial nearest neighbor kriging is used for the latent location-specific spatial parameter vectors \code{\alpha_{jl_j}}'s
#' to accelerate computation. Default is TRUE. When TRUE, \code{spatial.structure} must be \code{"continuous"}. If FALSE, the exact \code{M x M} neighborhood structure matrix \code{F(\rho)} 
#' is adopted throughout.
#'   
#' @param alphaMethod Character string indicating which approach is used to update the \code{\alpha_{jl_j}}'s in the Gibbs Sampler. 
#' Options include: \code{"block"} and \code{"sequential"}. This argument is only used when \code{include.space = TRUE} and \code{spatApprox = TRUE}.
#' 
#' @param h a positive integer much smaller than M, the number of location points, indicating an upper bound of the number of nearest neighbors for each location point. 
#' 
#' @param storeSpatPredPara a logical indicating whether we want to store the posterior samples for Xi, Delta, Tau, Alpha, Theta when clustering = TRUE.
#' Has to be TRUE if we want to perform predictions at new spatial location(s). Can be FALSE to save storage space if spatial prediction is not required.
#' 
#' @param storeWeights a logical indicating whether we want to store the posterior Weights estimates when clustering = TRUE. Has to be TRUE
#' if we want to perform temporal trends clustering from the output object. Can be FALSE to save storage space if otherwise.
#'
#' @param alphasWeightsToFiles a logical indicating whether we want to store the posterior alpha and weights estimates to files (one 
#' for each kept MCMC iteration) instead of including them in the main output list. Can only be TRUE if \code{clustering = TRUE} 
#' and at least one of \code{storeSpatPredPara} and \code{storeWeights} is TRUE. When that is the case, specifying \code{alphasWeightsToFiles = TRUE} 
#' would be immensely useful if we don't have enough space to output a model fit object with all kept posterior parameter estimates.
#' Default is FALSE. 
#' 
#'
#' @return \code{VAR1bfaFixedL} returns a list containing the following objects (some may be NULL)
#'
#'   \describe{
#'
#'   \item{\code{lambda}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for the factor loadings matrix \code{Lambda}.
#'   The labels for each column are Lambda_O_M_K.}
#'
#'   \item{\code{eta}}{\code{NKeep x (Nu x K)} \code{matrix} of posterior samples for the latent factors \code{eta}.
#'   The labels for each column are Eta_Nu_K.}
#'
#'   \item{\code{beta}}{\code{NKeep x P} \code{matrix} of posterior samples for \code{beta}.}
#'
#'   \item{\code{sigma2}}{\code{NKeep x (M * (O - C))} \code{matrix} of posterior samples for the variances \code{sigma2}.
#'   The labels for each column are Sigma2_O_M.}
#'
#'   \item{\code{kappa}}{\code{NKeep x ((O * (O + 1)) / 2)} \code{matrix} of posterior samples for \code{kappa}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Kappa3_2} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{delta}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{delta}.}
#'
#'   \item{\code{tau}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{tau}.}
#'
#'   \item{\code{upsilon}}{\code{NKeep x ((K * (K + 1)) / 2)} \code{matrix} of posterior samples for \code{Upsilon}. The
#'   columns have names that describe the samples within them. The row is listed first, e.g.,
#'   \code{Upsilon3_2} refers to the entry in row \code{3}, column \code{2}.}
#'
#'   \item{\code{A}}{\code{NKeep x (K * K)} \code{matrix} of posterior samples for \code{A}.}
#'
#'   \item{\code{xi}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for factor loadings cluster labels \code{xi}.
#'   The labels for each column are Xi_O_M_K.}
#'   
#'   \item{\code{rho}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{rho}.}
#'
#'   \item{\code{theta}}{\code{NKeep x (L x K)} \code{matrix} of posterior samples for \code{theta}.}
#'
#'   \item{\code{alpha}}{\code{NKeep x (M x O x K x (L - 1)} \code{matrix} of posterior samples for \code{alpha}.}
#'   For each kept row (MCMC iteration) in \code{alpha}, the corresponding entries are ordered first by observation type, then spatially, 
#'   then by clustering group, and finally by factor.
#'
#'   \item{\code{weights}}{\code{NKeep x (M x O x K x L)} \code{matrix} of posterior samples for \code{weights}.}
#'   For each kept row (MCMC iteration) in \code{weights}, the corresponding entries are ordered first spatially, then by observation type, 
#'   then by clustering group, and finally by factor.
#'
#'   \item{\code{metropolis}}{\code{2 (or 1) x 3} \code{matrix} of metropolis
#'   acceptance rates, updated tuners, and original tuners that result from the pilot
#'   adaptation.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{bfaFixedL} functions
#'   and should be ignored by the user.}
#'
#'   \item{\code{dataug}}{A \code{list} of data augmentation objects that are used in future
#'   \code{bfaFixedL} functions and should be ignored by the user.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'   
#'   }
#'
#' @export
VAR1bfaFixedL <- function(formula, data, dist, Nu, K, L = 20, 
                   family = "normal", spatial.structure = "continuous",
                   starting = NULL, hypers = NULL, tuning = NULL, mcmc = NULL, seed = 27,
                   gamma.shrinkage = TRUE, include.space = TRUE, clustering = TRUE, 
                   spatApprox = TRUE, alphaMethod = "block", h = 15,
                   storeSpatPredPara = TRUE, storeWeights = TRUE, alphasWeightsToFiles = FALSE) {

  ###Check for missing objects
  if (missing(formula)) stop("formula: missing")
  if (missing(data)) stop("data: missing")
  if (missing(dist)) stop("dist: missing")
  if (missing(Nu)) stop("T: missing")
  if (missing(K)) stop("K: missing")

  ###Check model inputs
  VAR1CheckInputsFixedL(formula, data, dist, Nu, K, L, starting, hypers, tuning, 
              mcmc, family, spatial.structure, gamma.shrinkage, 
              include.space, clustering, spatApprox, alphaMethod, h,
              storeSpatPredPara, storeWeights, alphasWeightsToFiles)

  ####Set seed for reproducibility
  set.seed(seed)

  ###Check if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- VAR1CreateDatObjFixedL(formula, data, dist, Nu, K, L, family, spatial.structure, 
                         gamma.shrinkage, include.space, clustering, spatApprox, alphaMethod, 
                         h, storeSpatPredPara, storeWeights, alphasWeightsToFiles)
  HyPara <- VAR1CreateHyPara(hypers, DatObj) 
  MetrObj <- VAR1CreateMetrObj(tuning, DatObj)
  Para <- VAR1CreateParaFixedL(starting, DatObj, HyPara)
  ParaCL <- CreateParaCLfixedL(starting, DatObj)
  Para <- CreateParaLambda(Para, DatObj, ParaCL)
  if (!(include.space & (spatial.structure == "continuous") & spatApprox)) {
    SpatPara <- CreateSpatPara1FixedL(starting, DatObj, HyPara)
  } else if (!clustering) {
    SpatPara <- CreateSpatPara2FixedL(starting, DatObj, HyPara)
  } else if (alphaMethod == "sequential"){
    SpatPara <- CreateSpatPara3(starting, DatObj, HyPara)
  } else if (alphaMethod == "block") {
    SpatPara <- CreateSpatPara2FixedL(starting, DatObj, HyPara)
  }
  DatAug <- CreateDatAug(DatObj)
  McmcObj <- CreateMcmc(mcmc, DatObj)
  RawSamples <- VAR1CreateStorageFixedL(DatObj, McmcObj)
  
  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run MCMC sampler in Rcpp
  RegObj <- VAR1bfaRcppFixedL(DatObj, HyPara, MetrObj, Para, ParaCL, SpatPara, DatAug, McmcObj, RawSamples, Interactive)

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Collect output to be returned
  if (include.space) {Metropolis <- VAR1SummarizeMetropolis(DatObj, MetrObj, RegObj$metropolis, McmcObj)} else {Metropolis <- NULL}
  Samples <- VAR1FormatSamplesFixedL(DatObj, RegObj$rawsamples)

  ###Return FixedLbfaVAR1 object
  FixedLbfaVAR1 <- list(lambda = Samples$Lambda,
                eta = Samples$Eta,
                beta = Samples$Beta,
                sigma2 = Samples$Sigma2,
                kappa = Samples$Kappa,
                delta = Samples$Delta,
                tau = Samples$Tau,
                upsilon = Samples$Upsilon,
                A = Samples$A,
                xi = Samples$Xi,
                rho = Samples$Rho,                
                theta = Samples$Theta,
                alpha = Samples$Alpha,
                weights = Samples$Weights,
                metropolis = Metropolis,
                datobj = DatObj,
                hypara = HyPara,
                dataug = DatAug,
                runtime = paste0("Model runtime: ",round(RunTime, 2), " ", attr(RunTime, "units")))
  FixedLbfaVAR1 <- structure(FixedLbfaVAR1, class = "FixedLbfaVAR1")
  return(FixedLbfaVAR1)

###End sampler
}
