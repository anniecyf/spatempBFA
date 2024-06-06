#' Spatial factor analysis using a Bayesian hierarchical model.
#'
#' \code{bfaVaryingLjs} is a Markov chain Monte Carlo (MCMC) sampler for a Bayesian spatial factor analysis model. The spatial component is 
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
#' @param time A \code{Nu} dimensional vector containing the observed time points in increasing order.
#'
#' @param K A scalar that indicates the dimension (i.e., quantity) of latent factors.
#'  
#' @param LjVec A vector of length \code{K} consisting of positive integers indicating the starting numbers of latent clusters for the K columns of the factor loadings matrix. 
#'  
#' @param starting Either \code{NULL} or a \code{list} containing starting values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the starting values may be specified.
#'
#'  When \code{NULL} is chosen then default starting values are automatically generated.
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon} or
#'  \code{Psi} containing appropriate objects. \code{Beta} (or \code{Delta}) must either be a \code{P} (or \code{K}) dimensional
#'  vector or a scalar (the scalar populates the entire vector). \code{Sigma2} must be either a \code{M x (O - C)} matrix or a scalar.
#'  \code{Kappa} must be a \code{O x O} dimensional matrix, \code{Rho} a scalar, \code{Upsilon} a \code{K x K} matrix, and \code{Psi} a scalar.
#'
#' @param hypers Either \code{NULL} or a \code{list} containing hyperparameter values
#'  to be specified for the MCMC sampler. If \code{NULL} is not chosen then none, some or all
#'  of the hyperparameter values may be specified.
#'
#'  When \code{NULL} is chosen then default hyperparameter values are automatically
#'  generated. These default hyperparameters are described in detail in (Berchuck et al.).
#'  Otherwise a \code{list} must be provided with names \code{Beta}, \code{Delta}, \code{Sigma2}, \code{Kappa}, \code{Rho}, \code{Upsilon} or
#'  \code{Psi} containing further hyperparameter information. These objects are themselves
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
#'  \code{Psi} is a \code{list} with two objects, dependent on if the temporal kernel is \code{exponential} or \code{ar1}.
#'  For \code{exponential}, the two objects are \code{APsi} and \code{BPsi}. \code{APsi}
#'  represents the lower bound for the uniform hyperprior, while \code{BPsi} represents
#'  the upper bound. The bounds must be specified carefully. For \code{ar1}, the two objects are \code{Beta} and \code{Gamma}, which are the 
#'  two shape parameters of a Beta distribution shifted to have domain in (-1, 1). 
#'  
#' @param tuning Either \code{NULL} or a \code{list} containing tuning values
#'  to be specified for the MCMC Metropolis steps. If \code{NULL} is not chosen then all
#'  of the tuning values must be specified.
#'
#'  When \code{NULL} is chosen then default tuning values are automatically generated to
#'  \code{1}. Otherwise a \code{list} must be provided with names \code{Psi}, 
#'  or \code{Rho}. Each of these entries must be scalars containing tuning variances for their corresponding Metropolis updates.
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
#' @param temporal.structure Character string indicating the temporal kernel. Options include:
#'  \code{"exponential"}, \code{"ar1"}, \code{"sar1"}, and \code{"sexponential"}.
#'
#' @param spatial.structure Character string indicating the type of spatial process. Options include:
#'  \code{"continuous"} (i.e., Gaussian process with exponential kernel) and \code{"discrete"} (i.e., proper CAR).
#'
#' @param seed An integer value used to set the seed for the random number generator
#'  
#' @param gamma.shrinkage A logical indicating whether a gamma shrinkage process prior is used for the variances of the factor loadings columns. If FALSE,
#'  the hyperparameters (A1 and A2) indicate the shape and rate for a gamma prior on the precisions. Default is TRUE. 
#'
#' @param include.time A logical indicating whether a temporal process should be included. Default is TRUE, however if FALSE the temporal correlation 
#'  structure matrix \code{H(\psi)} is fixed as the \code{T x T} identity matrix. This specification overrides all inputs \code{temporal.structure}, \code{equalTimeDist},
#'  and \code{seasonPeriod}.
#'
#' @param include.space A logical indicating whether a spatial process should be included. Default is TRUE, however if FALSE the spatial correlation matrix 
#'  \code{F(\rho)} is fixed as the \code{M x M} identity matrix. This specification overrides the inputs \code{spatial.structure}, \code{spatApprox}, and \code{alphaMethod}.
#'  
#' @param seasonPeriod A positive integer value indicating the temporal seasonality period. Default is 1 (no temporal seasonality). This argument is only used 
#'  (under which scenario it should be a positive integer greater than 1) when \code{include.time = TRUE} and the \code{temporal.structure} input corresponds 
#'  to a structure with seasonality.
#'
#' @param equalTimeDist A logical indicating whether the distances between adjacent time points are equal. Default is TRUE. If TRUE, the user should normalize 
#' the time distance to 1 (for the \code{time} input).
#'  
#' @param spatApprox A logical indicating whether spatial nearest neighbor kriging is used for the latent location-specific spatial parameter vectors \code{\alpha_{jl_j}}'s
#' to accelerate computation. Default is TRUE. When TRUE, \code{spatial.structure} must be \code{"continuous"}. If FALSE, the exact \code{M x M} neighborhood structure matrix \code{F(\rho)} 
#' is adopted throughout.
#' 
#' @param alphaSequen A logical indicating whether we want to sequentially update the alpha's. This argument is only used when \code{include.space = TRUE} and \code{spatApprox = TRUE}.
#' 
#' @param h a positive integer much smaller than M, the number of location points, indicating an upper bound of the number of nearest neighbors for each location point. 
#'
#' @param storeSpatPredPara a logical indicating whether we want to store the posterior samples for Alpha, Theta.
#' Has to be TRUE if we want to perform predictions at new spatial location(s). Can be FALSE to save storage space if spatial prediction is not required.
#' 
#' @param storeWeights a logical indicating whether we want to store the posterior Weights estimates. Has to be TRUE
#' if we want to perform temporal trends clustering from the output object. Can be FALSE to save storage space if otherwise.
#'
#' @return \code{bfaVaryingLjs} returns a list containing the following objects (some may be NULL)
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
#'   \item{\code{psi}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{psi}.}
#'
#'   \item{\code{xi}}{\code{NKeep x (M x O x K)} \code{matrix} of posterior samples for factor loadings cluster labels \code{xi}.
#'   The labels for each column are Xi_O_M_K.}
#'   
#'   \item{\code{rho}}{\code{NKeep x 1} \code{matrix} of posterior samples for \code{rho}.}
#'   
#'   \item{\code{ljvec}}{\code{NKeep x K} \code{matrix} of posterior samples for \code{Lj}'s.}
#'
#'   \item{\code{theta}}{A \code{list} of posterior samples for \code{theta}. The list contains \code{K x NKeep} elements (ordered by column), 
#'   each of which is a vector of length \code{Lj} at that MCMC iteration.}
#'
#'   \item{\code{alpha}}{A \code{list} of posterior samples for \code{alpha}. The list contains \code{K x NKeep} elements (ordered by column, i.e., the first K components of the list correspond to the first kept iteration and so on), 
#'   each of which is a \code{matrix} of dimension \code{(Lj - 1) x (M x O)} for \code{Lj} at that MCMC iteration. If \code{Lj = 1} for a certain \code{j}
#'   at a particular MCMC iteration, then that corresponding \code{matrix} in \code{alpha} is of dimension \code{1 x (M x O)} with entries all set to +Inf.
#'   In each element matrix for each clustering group, the corresponding entries are ordered first by observation type, then spatially.
#'   (the first O correspond to the first location point, the next O correspond to the second location point and so on)}
#'
#'   \item{\code{weights}}{A \code{list} of posterior samples for \code{weights}. The list contains \code{K x NKeep} elements (ordered by column), 
#'   each of which is a \code{matrix} of dimension \code{Lj x (M x O)} for \code{Lj} at that MCMC iteration.
#'   In each element matrix for each clustering group, the corresponding entries are ordered first spatially, then by observation type.
#'   (the first M correspond to the first observation type, the next M correspond to the second observation type and so on)}
#'
#'   \item{\code{metropolis}}{\code{2 (or 1) x 3} \code{matrix} of metropolis
#'   acceptance rates, updated tuners, and original tuners that result from the pilot
#'   adaptation.}
#'
#'   \item{\code{datobj}}{A \code{list} of data objects that are used in future \code{bfaVaryingLjs} functions
#'   and should be ignored by the user.}
#'
#'   \item{\code{dataug}}{A \code{list} of data augmentation objects that are used in future
#'   \code{bfaVaryingLjs} functions and should be ignored by the user.}
#'
#'   \item{\code{runtime}}{A \code{character} string giving the runtime of the MCMC sampler.}
#'   
#'   }
#'
#' @export
bfaVaryingLjs <- function(formula, data, dist, time, K, LjVec, 
                   family = "normal", temporal.structure = "exponential", spatial.structure = "continuous",
                   starting = NULL, hypers = NULL, tuning = NULL, mcmc = NULL, seed = 27,
                   gamma.shrinkage = TRUE, include.time = TRUE, include.space = TRUE, 
                   seasonPeriod = 1, equalTimeDist = TRUE, spatApprox = TRUE, alphaSequen = FALSE, h = 15,
                   storeSpatPredPara = TRUE, storeWeights = TRUE) {

  ###Check for missing objects
  if (missing(formula)) stop("formula: missing")
  if (missing(data)) stop("data: missing")
  if (missing(dist)) stop("dist: missing")
  if (missing(time)) stop("time: missing")
  if (missing(K)) stop("K: missing")

  ###Check model inputs
  CheckInputsVaryingLjs(formula, data, dist, time, K, LjVec, starting, hypers, tuning, 
              mcmc, family, temporal.structure, spatial.structure, 
              gamma.shrinkage, include.time, include.space,
              seasonPeriod, equalTimeDist, spatApprox, alphaSequen, h,
              storeSpatPredPara, storeWeights)

  ####Set seed for reproducibility
  set.seed(seed)

  ###Check if the job is interactive
  Interactive <- interactive()

  ###Create objects for use in sampler
  DatObj <- CreateDatObjVaryingLjs(formula, data, dist, time, K, family, 
                         temporal.structure, spatial.structure, 
                         gamma.shrinkage, include.time, include.space,
                         seasonPeriod, equalTimeDist, spatApprox, alphaSequen,
                         h, storeSpatPredPara, storeWeights)
  HyPara <- CreateHyPara(hypers, DatObj) 
  MetrObj <- CreateMetrObj(tuning, DatObj)
  Para <- CreateParaVaryingLjs(starting, DatObj, HyPara, LjVec)
  if (!(include.space & (spatial.structure == "continuous") & spatApprox)) {
    SpatPara <- CreateSpatPara1VaryingLjs(starting, DatObj, HyPara)
  } else if (!alphaSequen) {
    SpatPara <- CreateSpatPara2VaryingLjs(starting, DatObj, HyPara)
  } else {
    SpatPara <- CreateSpatPara3(starting, DatObj, HyPara)
  }
  DatAug <- CreateDatAug(DatObj)
  McmcObj <- CreateMcmc(mcmc, DatObj)
  RawSamples <- CreateStorageVaryingLjs(DatObj, McmcObj)
  
  ###Time MCMC sampler
  BeginTime <- Sys.time()

  ###Run MCMC sampler in Rcpp
  RegObj <- bfaRcppVaryingLjs(DatObj, HyPara, MetrObj, Para, SpatPara, DatAug, McmcObj, RawSamples, Interactive)
  if (storeSpatPredPara) {
    ThetaPost <- RegObj$thetafield
    AlphaPost <- RegObj$alphafield
  } else {
    ThetaPost <- NULL
    AlphaPost <- NULL
  }
  if (storeWeights) {
    WeightsPost <- RegObj$weightsfield
  } else {
    WeightsPost <- NULL
  }

  ###End time
  FinishTime <- Sys.time()
  RunTime <- FinishTime - BeginTime

  ###Collect output to be returned
  Metropolis <- SummarizeMetropolis(DatObj, MetrObj, RegObj$metropolis, McmcObj)
  Samples <- FormatSamplesVaryingLjs(DatObj, RegObj$rawsamples)

  ###Return varyLjBFA object
  varyLjBFA <- list(lambda = Samples$Lambda,
                eta = Samples$Eta,
                beta = Samples$Beta,
                sigma2 = Samples$Sigma2,
                kappa = Samples$Kappa,
                delta = Samples$Delta,
                tau = Samples$Tau,
                upsilon = Samples$Upsilon,
                psi = Samples$Psi,
                xi = Samples$Xi,
                rho = Samples$Rho,  
                ljvec = Samples$LjVec,
                theta = ThetaPost,
                alpha = AlphaPost,
                weights = WeightsPost,
                metropolis = Metropolis,
                datobj = DatObj,
                hypara = HyPara,
                dataug = DatAug,
                runtime = paste0("Model runtime: ",round(RunTime, 2), " ", attr(RunTime, "units")))
  varyLjBFA <- structure(varyLjBFA, class = "varyLjBFA")
  return(varyLjBFA)

###End sampler
}
