###Functions for reading in sampler inputs and creating a list object that contains all relevant data objects--------------------
CreateDatObjFixedL <- function(formula, data, dist, time, K, L, family, temporal.structure, spatial.structure, gamma.shrinkage, 
                         include.time, include.space, clustering, seasonPeriod, equalTimeDist, spatApprox, alphaMethod, h, 
                         storeSpatPredPara, storeWeights, alphasWeightsToFiles) {

  ###Data objects
  N <- dim(data)[1] # total observations
  M <- dim(dist)[1] #number of spatial locations
  Nu <- length(time) # number of visits
  O <- N / (M * Nu) #number of spatial observation types
  YObserved <- array(data[, all.vars(formula)[1]], dim = c(M, O, Nu)) # observed data
  
  if (!clustering) { L <- 1 }
  
  ###Covariates
  X <- model.matrix(formula, data = data)
  P <- dim(X)[2]
  X <- matrix(X, nrow = N, ncol = P)
  
  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(as.numeric(YObserved), ncol = 1)
  YStarWide <- matrix(YStar, nrow = M * O, ncol = Nu)
  
  ###Temporal distance matrix
  TimeDist <- abs(outer(time, time, "-"))

  ###Matrix Objects
  EyeNu <- diag(Nu)
  SeqL <- matrix(0:(L - 1), ncol = 1)
  Indeces <- matrix(rep(1:Nu - 1, each = M * O), ncol = 1)
  
  ###Assign temporal correlation structure
  if (temporal.structure == "exponential") TempCorInd <- 0
  if (temporal.structure == "ar1") TempCorInd <- 1
  if (temporal.structure == "sar1") TempCorInd <- 2
  if (temporal.structure == "sexponential") TempCorInd <- 3

  ###Assign spatial correlation structure
  if (spatial.structure == "continuous") {SpCorInd <- 0} else {SpCorInd <- 1} #if (spatial.structure == "discrete") 
  
  ###Approach to update the alpha's in the Gibbs Sampler when include.space = spatApprox = TRUE
  if (alphaMethod == "block") {alphaMethodInd <- 0} else {alphaMethodInd <- 1} #if (alphaMethod == "sequential") 
  
  ###Family indicator
  FamilyInd <- numeric(length = O)
  if (length(family) == 1) { # i.e., when either O = 1 or (O > 1 and length(family) = 1 all obs types same family)
    if (family == "normal") FamilyInd <- rep(0, O)
    if (family == "probit") FamilyInd <- rep(1, O)
    if (family == "tobit") FamilyInd <- rep(2, O)
    if (family == "binomial") FamilyInd <- rep(3, O)
  } else if (length(family) == O) { #i.e., when O > 1 and length(family) = O
    for (o in 1:O) {
      if (family[o] == "normal") FamilyInd[o] <- 0
      if (family[o] == "probit") FamilyInd[o] <- 1
      if (family[o] == "tobit") FamilyInd[o] <- 2
      if (family[o] == "binomial") FamilyInd[o] <- 3
    }
  }
    
  ###Assess the number of trials
  C <- sum(FamilyInd == 3)
  if (C > 0) {
    TrialsArray <- array(data$trials, dim = c(M, O, Nu))
    Trials <- array(dim = c(M, O, Nu))
    Trials[ , which(FamilyInd == 3), ] <- TrialsArray[ , which(FamilyInd == 3), ]
    Trials[ , which(FamilyInd != 3), ] <- 1
    Chi <- YObserved - 0.5 * Trials
  } else { #if (C == 0) 
    Chi <- Trials <- array(1, dim = c(1, 1, 1))
  }
  
  ###Make parameters global
  DatObj <- list()
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$SpDist <- dist
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$K <- K
  DatObj$L <- L
  DatObj$O <- O
  DatObj$C <- C
  DatObj$P <- P
  DatObj$X <- X
  DatObj$FamilyInd <- FamilyInd
  DatObj$Time <- time
  DatObj$TempCorInd <- TempCorInd
  DatObj$SpCorInd <- SpCorInd
  DatObj$alphaMethodInd <- alphaMethodInd
  DatObj$TimeDist <- TimeDist
  DatObj$EyeNu <- EyeNu
  DatObj$SeqL <- SeqL
  DatObj$Trials <- Trials
  DatObj$Chi <- Chi
  DatObj$Indeces <- Indeces
  DatObj$seasonPeriod <- seasonPeriod
  DatObj$h <- h
  DatObj$GS <- 1 * gamma.shrinkage
  DatObj$IT <- 1 * include.time
  DatObj$IS <- 1 * include.space
  DatObj$CL <- 1 * clustering
  DatObj$ET <- 1 * equalTimeDist
  DatObj$SA <- 1 * spatApprox
  DatObj$spatPred <- 1 * storeSpatPredPara
  DatObj$storeW <- 1 * storeWeights
  DatObj$alphaWsFiles <- 1 * alphasWeightsToFiles
  return(DatObj)

}

VAR1CreateDatObjFixedL <- function(formula, data, dist, Nu, K, L, family, spatial.structure, gamma.shrinkage, 
                               include.space, clustering, spatApprox, alphaMethod, h, 
                               storeSpatPredPara, storeWeights, alphasWeightsToFiles) {
  
  ###Data objects
  N <- dim(data)[1] # total observations
  M <- dim(dist)[1] #number of spatial locations
  O <- N / (M * Nu) #number of spatial observation types
  YObserved <- array(data[, all.vars(formula)[1]], dim = c(M, O, Nu)) # observed data
  
  if (!clustering) { L <- 1 }
  
  ###Covariates
  X <- model.matrix(formula, data = data)
  P <- dim(X)[2]
  X <- matrix(X, nrow = N, ncol = P)
  
  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(as.numeric(YObserved), ncol = 1)
  YStarWide <- matrix(YStar, nrow = M * O, ncol = Nu)

  ###Matrix Objects
  SeqL <- matrix(0:(L - 1), ncol = 1)
  Indeces <- matrix(rep(1:Nu - 1, each = M * O), ncol = 1)
  
  ###Assign spatial correlation structure
  if (spatial.structure == "continuous") {SpCorInd <- 0} else {SpCorInd <- 1} #if (spatial.structure == "discrete") 
  
  ###Approach to update the alpha's in the Gibbs Sampler when include.space = spatApprox = TRUE
  if (alphaMethod == "block") {alphaMethodInd <- 0} else {alphaMethodInd <- 1} #if (alphaMethod == "sequential")
  
  ###Family indicator
  FamilyInd <- numeric(length = O)
  if (length(family) == 1) { # i.e., when either O = 1 or (O > 1 and length(family) = 1 all obs types same family)
    if (family == "normal") FamilyInd <- rep(0, O)
    if (family == "probit") FamilyInd <- rep(1, O)
    if (family == "tobit") FamilyInd <- rep(2, O)
    if (family == "binomial") FamilyInd <- rep(3, O)
  } else if (length(family) == O) { #i.e., when O > 1 and length(family) = O
    for (o in 1:O) {
      if (family[o] == "normal") FamilyInd[o] <- 0
      if (family[o] == "probit") FamilyInd[o] <- 1
      if (family[o] == "tobit") FamilyInd[o] <- 2
      if (family[o] == "binomial") FamilyInd[o] <- 3
    }
  }
  
  ###Assess the number of trials
  C <- sum(FamilyInd == 3)
  if (C > 0) {
    TrialsArray <- array(data$trials, dim = c(M, O, Nu))
    Trials <- array(dim = c(M, O, Nu))
    Trials[ , which(FamilyInd == 3), ] <- TrialsArray[ , which(FamilyInd == 3), ]
    Trials[ , which(FamilyInd != 3), ] <- 1
    Chi <- YObserved - 0.5 * Trials
  } else { #if (C == 0)
    Chi <- Trials <- array(1, dim = c(1, 1, 1))
  }
  
  ###Make parameters global
  DatObj <- list()
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$SpDist <- dist
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$K <- K
  DatObj$L <- L
  DatObj$O <- O
  DatObj$C <- C
  DatObj$P <- P
  DatObj$X <- X
  DatObj$FamilyInd <- FamilyInd
  DatObj$SpCorInd <- SpCorInd
  DatObj$alphaMethodInd <- alphaMethodInd
  DatObj$SeqL <- SeqL
  DatObj$Trials <- Trials
  DatObj$Chi <- Chi
  DatObj$Indeces <- Indeces
  DatObj$h <- h
  DatObj$GS <- 1 * gamma.shrinkage
  DatObj$IS <- 1 * include.space
  DatObj$CL <- 1 * clustering
  DatObj$SA <- 1 * spatApprox
  DatObj$spatPred <- 1 * storeSpatPredPara
  DatObj$storeW <- 1 * storeWeights
  DatObj$alphaWsFiles <- 1 * alphasWeightsToFiles
  return(DatObj)
  
}


CreateDatObjVaryingLjs <- function(formula, data, dist, time, K, family, temporal.structure, spatial.structure, gamma.shrinkage, 
                         include.time, include.space, seasonPeriod, equalTimeDist, spatApprox, alphaSequen, h, storeSpatPredPara, storeWeights) {

  ###Data objects
  N <- dim(data)[1] # total observations
  M <- dim(dist)[1] #number of spatial locations
  Nu <- length(time) # number of visits
  O <- N / (M * Nu) #number of spatial observation types
  YObserved <- array(data[, all.vars(formula)[1]], dim = c(M, O, Nu)) # observed data
  
  ###Covariates
  X <- model.matrix(formula, data = data)
  P <- dim(X)[2]
  X <- matrix(X, nrow = N, ncol = P)
  
  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(as.numeric(YObserved), ncol = 1)
  YStarWide <- matrix(YStar, nrow = M * O, ncol = Nu)
  
  ###Temporal distance matrix
  TimeDist <- abs(outer(time, time, "-"))

  ###Matrix Objects
  EyeNu <- diag(Nu)
  Indeces <- matrix(rep(1:Nu - 1, each = M * O), ncol = 1)
  
  ###Assign temporal correlation structure
  if (temporal.structure == "exponential") TempCorInd <- 0
  if (temporal.structure == "ar1") TempCorInd <- 1
  if (temporal.structure == "sar1") TempCorInd <- 2
  if (temporal.structure == "sexponential") TempCorInd <- 3

  ###Assign spatial correlation structure
  if (spatial.structure == "continuous") {SpCorInd <- 0} else {SpCorInd <- 1} #if (spatial.structure == "discrete")
  
  ###Family indicator
  FamilyInd <- numeric(length = O)
  if (length(family) == 1) { # i.e., when either O = 1 or (O > 1 and length(family) = 1 all obs types same family)
    if (family == "normal") FamilyInd <- rep(0, O)
    if (family == "probit") FamilyInd <- rep(1, O)
    if (family == "tobit") FamilyInd <- rep(2, O)
    if (family == "binomial") FamilyInd <- rep(3, O)
  } else if (length(family) == O) { #i.e., when O > 1 and length(family) = O
    for (o in 1:O) {
      if (family[o] == "normal") FamilyInd[o] <- 0
      if (family[o] == "probit") FamilyInd[o] <- 1
      if (family[o] == "tobit") FamilyInd[o] <- 2
      if (family[o] == "binomial") FamilyInd[o] <- 3
    }
  }
    
  ###Assess the number of trials
  C <- sum(FamilyInd == 3)
  if (C > 0) {
    TrialsArray <- array(data$trials, dim = c(M, O, Nu))
    Trials <- array(dim = c(M, O, Nu))
    Trials[ , which(FamilyInd == 3), ] <- TrialsArray[ , which(FamilyInd == 3), ]
    Trials[ , which(FamilyInd != 3), ] <- 1
    Chi <- YObserved - 0.5 * Trials
  } else { #if (C == 0) 
    Chi <- Trials <- array(1, dim = c(1, 1, 1)) 
  }

  ###Make parameters global
  DatObj <- list()
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$SpDist <- dist
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$K <- K
  DatObj$O <- O
  DatObj$C <- C
  DatObj$P <- P
  DatObj$X <- X
  DatObj$FamilyInd <- FamilyInd
  DatObj$Time <- time
  DatObj$TempCorInd <- TempCorInd
  DatObj$SpCorInd <- SpCorInd
  DatObj$TimeDist <- TimeDist
  DatObj$EyeNu <- EyeNu
  DatObj$Trials <- Trials
  DatObj$Chi <- Chi
  DatObj$Indeces <- Indeces
  DatObj$seasonPeriod <- seasonPeriod
  DatObj$h <- h
  DatObj$GS <- 1 * gamma.shrinkage
  DatObj$IT <- 1 * include.time
  DatObj$IS <- 1 * include.space
  DatObj$ET <- 1 * equalTimeDist
  DatObj$SA <- 1 * spatApprox
  DatObj$alphaSequenInd <- 1 * alphaSequen
  DatObj$spatPred <- 1 * storeSpatPredPara
  DatObj$storeW <- 1 * storeWeights
  return(DatObj)

}

VAR1CreateDatObjVaryingLjs <- function(formula, data, dist, Nu, K, family, spatial.structure, gamma.shrinkage, 
                                   include.space, spatApprox, alphaSequen, h, storeSpatPredPara, storeWeights) {
  
  ###Data objects
  N <- dim(data)[1] # total observations
  M <- dim(dist)[1] #number of spatial locations
  O <- N / (M * Nu) #number of spatial observation types
  YObserved <- array(data[, all.vars(formula)[1]], dim = c(M, O, Nu)) # observed data
  
  ###Covariates
  X <- model.matrix(formula, data = data)
  P <- dim(X)[2]
  X <- matrix(X, nrow = N, ncol = P)
  
  ###Dynamic Objects (updated with data augmentation)
  YStar <- matrix(as.numeric(YObserved), ncol = 1)
  YStarWide <- matrix(YStar, nrow = M * O, ncol = Nu)
  
  ###Matrix Objects
  Indeces <- matrix(rep(1:Nu - 1, each = M * O), ncol = 1)
  
  ###Assign spatial correlation structure
  if (spatial.structure == "continuous") {SpCorInd <- 0} else {SpCorInd <- 1} #if (spatial.structure == "discrete") 
  
  ###Family indicator
  FamilyInd <- numeric(length = O)
  if (length(family) == 1) { # i.e., when either O = 1 or (O > 1 and length(family) = 1 all obs types same family)
    if (family == "normal") FamilyInd <- rep(0, O)
    if (family == "probit") FamilyInd <- rep(1, O)
    if (family == "tobit") FamilyInd <- rep(2, O)
    if (family == "binomial") FamilyInd <- rep(3, O)
  } else if (length(family) == O) { #i.e., when O > 1 and length(family) = O
    for (o in 1:O) {
      if (family[o] == "normal") FamilyInd[o] <- 0
      if (family[o] == "probit") FamilyInd[o] <- 1
      if (family[o] == "tobit") FamilyInd[o] <- 2
      if (family[o] == "binomial") FamilyInd[o] <- 3
    }
  }
  
  ###Assess the number of trials
  C <- sum(FamilyInd == 3)
  if (C > 0) {
    TrialsArray <- array(data$trials, dim = c(M, O, Nu))
    Trials <- array(dim = c(M, O, Nu))
    Trials[ , which(FamilyInd == 3), ] <- TrialsArray[ , which(FamilyInd == 3), ]
    Trials[ , which(FamilyInd != 3), ] <- 1
    Chi <- YObserved - 0.5 * Trials
  } else {  #if (C == 0) 
    Chi <- Trials <- array(1, dim = c(1, 1, 1))
  }
  
  ###Make parameters global
  DatObj <- list()
  DatObj$YStar <- YStar
  DatObj$YStarWide <- YStarWide
  DatObj$SpDist <- dist
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$K <- K
  DatObj$O <- O
  DatObj$C <- C
  DatObj$P <- P
  DatObj$X <- X
  DatObj$FamilyInd <- FamilyInd
  DatObj$SpCorInd <- SpCorInd
  DatObj$Trials <- Trials
  DatObj$Chi <- Chi
  DatObj$Indeces <- Indeces
  DatObj$h <- h
  DatObj$GS <- 1 * gamma.shrinkage
  DatObj$IS <- 1 * include.space
  DatObj$SA <- 1 * spatApprox
  DatObj$alphaSequenInd <- 1 * alphaSequen
  DatObj$spatPred <- 1 * storeSpatPredPara
  DatObj$storeW <- 1 * storeWeights
  return(DatObj)
  
}



###Function to create Hyperparameter Object------------------------------------------------------------------------------------
CreateHyPara <- function(hypers, DatObj) {

  ###Set data objects
  K <- DatObj$K
  O <- DatObj$O
  P <- DatObj$P
  TempCorInd <- DatObj$TempCorInd
  SpCorInd <- DatObj$SpCorInd
  TimeDist <- DatObj$TimeDist 
  SpDist <- DatObj$SpDist

  ###Which parameters are user defined?
  Userhypers <- names(hypers)

  ###Set hyperparameters for Sigma2
  if ("Sigma2" %in% Userhypers) {
    A <- hypers$Sigma2$A
    B <- hypers$Sigma2$B
  } else {
    A <- 1
    B <- 1
  }
  
  ###Set hyperparameters for Kappa
  if ("Kappa" %in% Userhypers) {
    SmallUpsilon <- hypers$Kappa$SmallUpsilon
    BigTheta <- matrix(hypers$Kappa$BigTheta, nrow = O, ncol = O)
  } else {
    if (O > 1) {SmallUpsilon <- O + 1} else {SmallUpsilon <- 0.001}#if O == 1
    if (O == 1) {
        BigTheta <- 0.001 * diag(O)
    } else if (O <= 4) {
        BigTheta <- 0.1 * diag(O)
    } else {
        BigTheta <- diag(O)
    }
  }

  ###Set hyperparameters for Rho
  if ("Rho" %in% Userhypers) {
    if (SpCorInd == 0) { # continuous
      ARho <- hypers$Rho$ARho
      BRho <- hypers$Rho$BRho
    } else { # discrete (if (SpCorInd == 1) )
      ARho <- 0 
      BRho <- 1
    }
  } else {
    if (SpCorInd == 0) { # continuous
      minDiff <- min(SpDist[SpDist > 0])
      maxDiff <- max(SpDist[SpDist > 0])
      ARho <- -log(0.95) / maxDiff #longest diff goes up to 95%
      BRho <- -log(0.01) / minDiff #shortest diff goes down to 1%
    } else { # discrete (if (SpCorInd == 1) )
      ARho <- 0 
      BRho <- 1
    }
  }
  
  ###Set hyperparameters for Delta
  if ("Delta" %in% Userhypers) {
    A1 <- hypers$Delta$A1
    A2 <- hypers$Delta$A2
  } else {
    A1 <- 1
    A2 <- 1
  }

  ###Set hyperparameters for Beta
  if ("Beta" %in% Userhypers) {
    MuBeta <- hypers$Beta$MuBeta
    SigmaBeta <- hypers$Beta$SigmaBeta
    SigmaBetaInv <- CholInv(SigmaBeta)
    SigmaBetaInvMuBeta <- SigmaBetaInv %*% MuBeta
  } else {
    MuBeta <- matrix(0, nrow = P, ncol = 1)
    SigmaBetaInv <- diag(rep(1 / 1000, P))
    SigmaBetaInvMuBeta <- SigmaBetaInv %*% MuBeta
  }
  
  ###Set hyperparameters for Psi
  if ("Psi" %in% Userhypers) {
    if (TempCorInd == 0) { # exponential
      APsi <- hypers$Psi$APsi
      BPsi <- hypers$Psi$BPsi
      Gamma <- 1 #Null values
      Beta <- 1 #Null values
    } else if (TempCorInd == 1) { # ar1
      APsi <- -1
      BPsi <- 1
      Gamma <- hypers$Psi$Gamma
      Beta <- hypers$Psi$Beta
    } else if (TempCorInd == 2){ # sar1
      APsi <- -1
      BPsi <- 1
      Gamma <- hypers$Psi$Gamma
      Beta <- hypers$Psi$Beta
    } else if (TempCorInd == 3){ # sexponential
      APsi <- hypers$Psi$APsi
      BPsi <- hypers$Psi$BPsi
      Gamma <- 1 #Null values
      Beta <- 1 #Null values
    }
  } else { #if (!"Psi" %in% Userhypers) 
    if (TempCorInd == 0) { # exponential
      minDiff <- min(TimeDist[TimeDist > 0])
      maxDiff <- max(TimeDist[TimeDist > 0])
      APsi <- -log(0.95) / maxDiff #longest diff goes up to 95%
      BPsi <- -log(0.01) / minDiff #shortest diff goes down to 1%
      Gamma <- 1 #Null values
      Beta <- 1 #Null values
    } else if (TempCorInd == 1) { # ar1
      APsi <- -1
      BPsi <- 1
      Gamma <- 1
      Beta <- 1
    } else if (TempCorInd == 2) { # sar1
      APsi <- -1
      BPsi <- 1
      Gamma <- 1
      Beta <- 1
    } else if (TempCorInd == 3) { # sexponential
      minDiff <- min(TimeDist[TimeDist > 0])
      maxDiff <- max(TimeDist[TimeDist > 0])
      APsi <- -log(0.95) / maxDiff #longest diff goes up to 95%
      BPsi <- -log(0.01) / minDiff #shortest diff goes down to 1%
      Gamma <- 1 #Null values
      Beta <- 1 #Null values
    }
  }
  
  ###Set hyperparameters for Upsilon
  if ("Upsilon" %in% Userhypers) {
    Zeta <- hypers$Upsilon$Zeta
    Omega <- matrix(hypers$Upsilon$Omega, nrow = K, ncol = K)
  } else { # if (!("Upsilon" %in% Userhypers)) 
    if (K > 1) {Zeta <- K + 1} else {Zeta <- 0.001}
    if (K == 1) {
        Omega <- 0.001 * diag(K)
    } else if (K <= 4) {
        Omega <- 0.1 * diag(K)
    } else {
        Omega <- diag(K)
    }
  }

  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$A <- A
  HyPara$B <- B
  HyPara$SmallUpsilon <- SmallUpsilon
  HyPara$BigTheta <- BigTheta
  HyPara$A1 <- A1
  HyPara$A2 <- A2
  HyPara$APsi <- APsi
  HyPara$BPsi <- BPsi
  HyPara$Gamma <- Gamma
  HyPara$Beta <- Beta
  HyPara$Zeta <- Zeta
  HyPara$Omega <- Omega
  HyPara$ARho <- ARho
  HyPara$BRho <- BRho
  HyPara$SigmaBetaInvMuBeta <- SigmaBetaInvMuBeta
  HyPara$SigmaBetaInv <- SigmaBetaInv
  return(HyPara)

}

VAR1CreateHyPara <- function(hypers, DatObj) {
  
  ###Set data objects
  K <- DatObj$K
  O <- DatObj$O
  P <- DatObj$P
  SpCorInd <- DatObj$SpCorInd
  SpDist <- DatObj$SpDist
  
  ###Which parameters are user defined?
  Userhypers <- names(hypers)
  
  ###Set hyperparameters for Sigma2
  if ("Sigma2" %in% Userhypers) {
    A <- hypers$Sigma2$A
    B <- hypers$Sigma2$B
  } else {
    A <- 1
    B <- 1
  }
  
  ###Set hyperparameters for Kappa
  if ("Kappa" %in% Userhypers) {
    SmallUpsilon <- hypers$Kappa$SmallUpsilon
    BigTheta <- matrix(hypers$Kappa$BigTheta, nrow = O, ncol = O)
  } else {
    if (O > 1) {SmallUpsilon <- O + 1} else {SmallUpsilon <- 0.001}#if O == 1
    if (O == 1) {
      BigTheta <- 0.001 * diag(O)
    } else if (O <= 4){
      BigTheta <- 0.1 * diag(O)
    } else {
      BigTheta <- diag(O)
    }
  }
  
  ###Set hyperparameters for Rho
  if ("Rho" %in% Userhypers) {
    if (SpCorInd == 0) { # continuous
      ARho <- hypers$Rho$ARho
      BRho <- hypers$Rho$BRho
    } else { # discrete (if (SpCorInd == 1) )
      ARho <- 0 
      BRho <- 1
    }
  } else {
    if (SpCorInd == 0) { # continuous
      minDiff <- min(SpDist[SpDist > 0])
      maxDiff <- max(SpDist[SpDist > 0])
      ARho <- -log(0.95) / maxDiff #longest diff goes up to 95%
      BRho <- -log(0.01) / minDiff #shortest diff goes down to 1%
    } else { # discrete (if (SpCorInd == 1) )
      ARho <- 0 
      BRho <- 1
    }
  }
  
  ###Set hyperparameters for Delta
  if ("Delta" %in% Userhypers) {
    A1 <- hypers$Delta$A1
    A2 <- hypers$Delta$A2
  } else {
    A1 <- 1
    A2 <- 1
  }
  
  ###Set hyperparameters for Beta
  if ("Beta" %in% Userhypers) {
    MuBeta <- hypers$Beta$MuBeta
    SigmaBeta <- hypers$Beta$SigmaBeta
    SigmaBetaInv <- CholInv(SigmaBeta)
    SigmaBetaInvMuBeta <- SigmaBetaInv %*% MuBeta
  } else {
    MuBeta <- matrix(0, nrow = P, ncol = 1)
    SigmaBetaInv <- diag(rep(1 / 1000, P))
    SigmaBetaInvMuBeta <- SigmaBetaInv %*% MuBeta
  }
  
  ###Set hyperparameters for Upsilon
  if ("Upsilon" %in% Userhypers) {
    Zeta <- hypers$Upsilon$Zeta
    Omega <- matrix(hypers$Upsilon$Omega, nrow = K, ncol = K)
  } else { # if (!("Upsilon" %in% Userhypers)) 
    if (K > 1) {Zeta <- K + 1} else {Zeta <- 0.001}
    if (K == 1) {
      Omega <- 0.001 * diag(K)
    } else if (K <= 4) {
      Omega <- 0.1 * diag(K)
    } else {
      Omega <- diag(K)
    }
  }
  
  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$A <- A
  HyPara$B <- B
  HyPara$SmallUpsilon <- SmallUpsilon
  HyPara$BigTheta <- BigTheta
  HyPara$A1 <- A1
  HyPara$A2 <- A2
  HyPara$Zeta <- Zeta
  HyPara$Omega <- Omega
  HyPara$ARho <- ARho
  HyPara$BRho <- BRho
  HyPara$SigmaBetaInvMuBeta <- SigmaBetaInvMuBeta
  HyPara$SigmaBetaInv <- SigmaBetaInv
  return(HyPara)
  
}



###Function for creating an object containing relevant Metropolis information---------------------------------------------------
CreateMetrObj <- function(tuning, DatObj) {
  
  ###Which parameters are user defined?
  UserTuners <- names(tuning)

  ###Set tuning parameters for Psi and Rho
  if ("Psi" %in% UserTuners) {MetropPsi <- tuning$Psi} else {MetropPsi <- 1}
  if ("Rho" %in% UserTuners) { MetropRho <- tuning$Rho } else { MetropRho <- 1 }
  
  ###Set acceptance rate counters
  AcceptancePsi <- AcceptanceRho <- 0

  ###Return metropolis object
  MetrObj <- list()
  MetrObj$MetropPsi <- MetropPsi
  MetrObj$AcceptancePsi <- AcceptancePsi
  MetrObj$MetropRho <- MetropRho
  MetrObj$AcceptanceRho <- AcceptanceRho
  MetrObj$OriginalTuners <- c(MetropPsi, MetropRho)
  return(MetrObj)

}

VAR1CreateMetrObj <- function(tuning, DatObj) {
  
  UserTuners <- names(tuning)
  if ("Rho" %in% UserTuners) { MetropRho <- tuning$Rho } else { MetropRho <- 1 }
  
  AcceptanceRho <- 0
  MetrObj <- list()
  MetrObj$MetropRho <- MetropRho
  MetrObj$AcceptanceRho <- AcceptanceRho
  MetrObj$OriginalTuner <- MetropRho
  return(MetrObj)
  
}



###Functions for creating initial parameters 
## corresponding to the fixedL main function (except for spatial ones and clustering ones)
CreateParaFixedL <- function(starting, DatObj, HyPara) {

  ###Set data objects
  K <- DatObj$K
  L <- DatObj$L
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  TempCorInd <- DatObj$TempCorInd
  TimeDist <- DatObj$TimeDist
  EyeNu <- DatObj$EyeNu
  FamilyInd <- DatObj$FamilyInd
  X <- DatObj$X
  IT <- DatObj$IT
  seasonPeriod <- DatObj$seasonPeriod

  ###Set hyperparameter objects
  APsi <- HyPara$APsi
  BPsi <- HyPara$BPsi
  
  ###Which parameters are user defined?
  UserStarters <- names(starting)

  ###Set initial values of Sigma2
  if (any(FamilyInd != 3)) {
    if ("Sigma2" %in% UserStarters) {Sigma2 <- matrix(starting$Sigma2, ncol = (O - C), nrow = M)} else {Sigma2 <- matrix(1, ncol = (O - C), nrow = M)}
  } else {Sigma2 <- matrix(1, nrow = 1, ncol = 1)}
  
  ###Set initial values of Kappa
  if ("Kappa" %in% UserStarters) {Kappa <- matrix(starting$Kappa, nrow = O, ncol = O)} else {Kappa <- diag(O)}
  
  ###Set initial values of Psi and Temporal parameters
  if (IT == 0){
    HPsi <- EyeNu
    RootiHPsi <- EyeNu
    HPsiInv <- EyeNu
    Psi <- 0
  } else { #if (IT == 1)
    if ("Psi" %in% UserStarters) {
      Psi <- starting$Psi
      if ((Psi <= APsi) | (Psi >= BPsi)) stop('starting: "Psi" must be in (APsi, BPsi)')
    } else { #if ((!"Psi" %in% UserStarters)) 
      if (TempCorInd == 0) { #exponential
        Psi <- mean(c(APsi, BPsi))
      } else if (TempCorInd == 1) { #ar1
        Psi <- 0.7
      } else if (TempCorInd == 2){ #sar1
        Psi <- 0.7
      } else if (TempCorInd == 3){ #sexponential
        Psi <- mean(c(APsi, BPsi))
      }
    }
    ###Temporal parameters
    HPsi <- getH(Psi, TempCorInd, TimeDist, Nu, seasonPeriod)
    CholHPsi <- chol(HPsi)
    RootiHPsi <- solve(CholHPsi)
    HPsiInv <- RootiHPsi %*% t(RootiHPsi)
    #HPsiInv <- chol2inv(CholHPsi)
  }
  
  
  ###Set initial values of Beta
  if ("Beta" %in% UserStarters) {Beta <- matrix(starting$Beta, nrow = P, ncol = 1)} else {Beta <- matrix(0, nrow = P, ncol = 1)}
    
  ###Set initial value of Upsilon
  if ("Upsilon" %in% UserStarters) {Upsilon <- matrix(starting$Upsilon, nrow = K, ncol = K)} else {Upsilon <- diag(K)}
  RootiUpsilon = solve(chol(Upsilon))
  UpsilonInv = RootiUpsilon %*% t(RootiUpsilon)

  ###Factors
  BigPhi <- matrix(0, nrow = K, ncol = Nu)
  Eta <- matrix(as.numeric(BigPhi), ncol = 1)
  
  ###Spatial covariance objects
  CholKappa <- chol(Kappa)
  RootiKappa <- solve(CholKappa)
  KappaInv <- RootiKappa %*% t(RootiKappa)
  #KappaInv <- chol2inv(CholKappa)

  ###Create Poly-Gamma update objects
  Cov <- array(dim = c(M, O, Nu))
  NotCount <- which(FamilyInd != 3)
  count <- 1
  for (f in 1:O) {
    if (FamilyInd[f] %in% 0:2) {
      Cov[ , f, ] <- matrix(Sigma2[, count], nrow = M, ncol = Nu)
      count <- count + 1
    }
    if (FamilyInd[f] == 3) {
      omega <- matrix(pgdraw(rep(1, M*Nu), rep(0, M*Nu)), nrow = M, ncol = Nu)
      Cov[ , f, ] <- 1 / omega
    }
  }

  ###Save parameter objects
  Para <- list()
  Para$Sigma2 <- Sigma2
  Para$Kappa <- Kappa
  Para$Beta <- Beta
  Para$XBeta <- X %*% Beta
  Para$Psi <- Psi
  Para$Upsilon <- Upsilon
  Para$UpsilonInv <- UpsilonInv
  Para$RootiUpsilon <- RootiUpsilon
  Para$BigPhi <- BigPhi
  Para$Eta <- Eta
  Para$HPsi <- HPsi
  Para$RootiHPsi <- RootiHPsi
  Para$HPsiInv <- HPsiInv
  Para$RootiKappa <- RootiKappa
  Para$KappaInv <- KappaInv
  Para$Cov <- Cov
  return(Para)
}

VAR1CreateParaFixedL <- function(starting, DatObj, HyPara) {
  
  ###Set data objects
  K <- DatObj$K
  L <- DatObj$L
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  FamilyInd <- DatObj$FamilyInd
  X <- DatObj$X
  
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  ###Set initial values of Sigma2
  if (any(FamilyInd != 3)) {
    if ("Sigma2" %in% UserStarters) {Sigma2 <- matrix(starting$Sigma2, ncol = (O - C), nrow = M)} else {Sigma2 <- matrix(1, ncol = (O - C), nrow = M)}
  } else {Sigma2 <- matrix(1, nrow = 1, ncol = 1)}
  
  ###Set initial values of Kappa
  if ("Kappa" %in% UserStarters) {Kappa <- matrix(starting$Kappa, nrow = O, ncol = O)} else {Kappa <- diag(O)}
  
  ###Set initial values of A
  if ("A" %in% UserStarters) { A <- starting$A } else { A <- 0.8 * diag(K) }
  
  ###Set initial values of Beta
  if ("Beta" %in% UserStarters) {Beta <- matrix(starting$Beta, nrow = P, ncol = 1)} else {Beta <- matrix(0, nrow = P, ncol = 1)}
  
  ###Set initial value of Upsilon
  if ("Upsilon" %in% UserStarters) {Upsilon <- matrix(starting$Upsilon, nrow = K, ncol = K)} else {Upsilon <- diag(K)}
  RootiUpsilon = solve(chol(Upsilon))
  UpsilonInv = RootiUpsilon %*% t(RootiUpsilon)
  
  ###Factors
  BigPhi <- matrix(0, nrow = K, ncol = Nu)
  Eta <- matrix(as.numeric(BigPhi), ncol = 1)
  
  ###Spatial covariance objects
  CholKappa <- chol(Kappa)
  RootiKappa <- solve(CholKappa)
  KappaInv <- RootiKappa %*% t(RootiKappa)
  #KappaInv <- chol2inv(CholKappa)
  
  ###Create Poly-Gamma update objects
  Cov <- array(dim = c(M, O, Nu))
  NotCount <- which(FamilyInd != 3)
  count <- 1
  for (f in 1:O) {
    if (FamilyInd[f] %in% 0:2) {
      Cov[ , f, ] <- matrix(Sigma2[, count], nrow = M, ncol = Nu)
      count <- count + 1
    }
    if (FamilyInd[f] == 3) {
      omega <- matrix(pgdraw(rep(1, M*Nu), rep(0, M*Nu)), nrow = M, ncol = Nu)
      Cov[ , f, ] <- 1 / omega
    }
  }
  
  ###Save parameter objects
  Para <- list()
  Para$Sigma2 <- Sigma2
  Para$Kappa <- Kappa
  Para$Beta <- Beta
  Para$XBeta <- X %*% Beta
  Para$A <- A
  Para$Upsilon <- Upsilon
  Para$UpsilonInv <- UpsilonInv
  Para$BigPhi <- BigPhi
  Para$Eta <- Eta
  Para$RootiKappa <- RootiKappa
  Para$KappaInv <- KappaInv
  Para$Cov <- Cov
  return(Para)
}

## corresponding to the fixedL main function when CL = 1 (otherwise returns an empty list)
CreateParaCLfixedL <- function(starting, DatObj){
  ParaCL <- list()
  
  if (DatObj$CL == 1){
    K <- DatObj$K
    L <- DatObj$L
    M <- DatObj$M
    O <- DatObj$O
    GS <- DatObj$GS
    ###Set initial values of Delta
    if ("Delta" %in% names(starting)) {Delta <- matrix(starting$Delta, nrow = K, ncol = 1)} else {Delta <- matrix(1, nrow = K, ncol = 1)} #if ((!"Delta" %in% names(starting))) 
    ###Create atom variances (Tau2) and atoms themselves (Theta)
    if (GS == 1) {Tau <- matrix(cumprod(Delta), nrow = K, ncol = 1)} else {Tau <- Delta} #if (GS == 0) 
    Theta <- matrix(sapply(sqrt(1/Tau), rnorm, n=L, mean=0), nrow=L, K) #L by K matrix
    ###Create label parameters
    Xi <- matrix(0, nrow = M * O, ncol = K)
    ###Probit parameters
    Alpha <- array(0, dim = c(L - 1, M * O, K))
    Weights <- logWeights <- array(0, dim = c(L, M * O, K))
    Z <- array(-1, dim = c(L - 1, M * O, K))
    Z[1, , ] <- 1
    logUpperPhiAlpha <- pnorm(Alpha, lower.tail = FALSE, log.p = TRUE)
    UpperPhiAlpha <- pnorm(Alpha, lower.tail = FALSE)
    for (j in 1:K) {
      for (o in 1:O) {
        for (i in 1:M) {
          Index1 <- i + M * (o - 1)
          Index2 <- o + O * (i - 1)
          for (l in 1:L) {
            if (l == 1) {
              Weights[l, Index1, j] <- pnorm(Alpha[l, Index2, j])
              logWeights[l, Index1, j] <- pnorm(Alpha[l, Index2, j], log.p = TRUE)
            } else if (l < L) {
              Weights[l, Index1, j] <- pnorm(Alpha[l, Index2, j]) * prod(UpperPhiAlpha[1:(l - 1) , Index2, j])
              logWeights[l, Index1, j] <- pnorm(Alpha[l, Index2, j], log.p = TRUE) + sum(logUpperPhiAlpha[1:(l - 1), Index2, j])
            } else {
              Weights[l, Index1, j] <- prod(UpperPhiAlpha[1:(l - 1) , Index2, j])
              logWeights[l, Index1, j] <- sum(logUpperPhiAlpha[1:(l - 1), Index2, j])
            }
          }
        } 
      }
    }
    if (DatObj$storeW == 0) {Weights <- array(1, dim = c(1, 1, 1))}
    ParaCL$Delta <- Delta
    ParaCL$Tau <- Tau
    ParaCL$Xi <- Xi
    ParaCL$Theta <- Theta
    ParaCL$Alpha <- Alpha
    ParaCL$Z <- Z
    ParaCL$Weights <- Weights
    ParaCL$logWeights <- logWeights
  }
  
  return(ParaCL)
}

## Augment Lambda to the Para list
CreateParaLambda <- function(Para, DatObj, ParaCL){
  K <- DatObj$K
  M <- DatObj$M
  O <- DatObj$O
  Lambda <- matrix(0, nrow = M * O, ncol = K)
  XBeta <- Para$XBeta
  Nu <- DatObj$Nu
  EyeNu <- diag(Nu)
  Eta <- Para$Eta
  
  if (DatObj$CL == 1){
    Theta <- ParaCL$Theta
    Xi <- ParaCL$Xi
    for (o in 1:O) {
      for (i in 1:M) {
        for (j in 1:K) {
          Lambda[i + M * (o - 1), j] <- Theta[Xi[i + M * (o - 1), j] + 1, j]
        }
      }  
    }
  }
  
  Para$Lambda <- Lambda
  Para$Mean <- kronecker(EyeNu, Lambda) %*% Eta + XBeta
  return(Para)
}


## corresponding to the VaryingLjs main function (for parameters with fixed dimensions only, i.e., excluding Theta, Alpha, Weights) 
CreateParaVaryingLjs <- function(starting, DatObj, HyPara, LjVec) {

  ###Set data objects
  K <- DatObj$K
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  TempCorInd <- DatObj$TempCorInd
  TimeDist <- DatObj$TimeDist
  EyeNu <- DatObj$EyeNu
  FamilyInd <- DatObj$FamilyInd
  X <- DatObj$X
  GS <- DatObj$GS
  IT <- DatObj$IT
  seasonPeriod <- DatObj$seasonPeriod

  ###Set hyperparameter objects
  APsi <- HyPara$APsi
  BPsi <- HyPara$BPsi
  
  ###Which parameters are user defined?
  UserStarters <- names(starting)

  ###Set initial values of Sigma2
  if (any(FamilyInd != 3)) {
    if ("Sigma2" %in% UserStarters) {Sigma2 <- matrix(starting$Sigma2, ncol = (O - C), nrow = M)} else {Sigma2 <- matrix(1, ncol = (O - C), nrow = M)} #if ((!"Sigma2" %in% UserStarters)) 
  } else {Sigma2 <- matrix(1, nrow = 1, ncol = 1)}
  
  ###Set initial values of Kappa
  if ("Kappa" %in% UserStarters) {Kappa <- matrix(starting$Kappa, nrow = O, ncol = O)} else {Kappa <- diag(O)} #if ((!"Kappa" %in% UserStarters)) 
  
  ###Set initial values of Psi and Temporal parameters
  if (IT == 0){
    HPsi <- EyeNu
    RootiHPsi <- EyeNu
    HPsiInv <- EyeNu
    Psi <- 0
  } else { #if (IT == 1)
    if ("Psi" %in% UserStarters) {
      Psi <- starting$Psi
      if ((Psi <= APsi) | (Psi >= BPsi)) stop('starting: "Psi" must be in (APsi, BPsi)')
    } else if ((!"Psi" %in% UserStarters)) {
      if (TempCorInd == 0) { #exponential
        Psi <- mean(c(APsi, BPsi))
      } else if (TempCorInd == 1) { #ar1
        Psi <- 0.7
      } else if (TempCorInd == 2){ #sar1
        Psi <- 0.7
      } else if (TempCorInd == 3){ #sexponential
        Psi <- mean(c(APsi, BPsi))
      }
    }
    ###Temporal parameters
    HPsi <- getH(Psi, TempCorInd, TimeDist, Nu, seasonPeriod)
    CholHPsi <- chol(HPsi)
    RootiHPsi <- solve(CholHPsi)
    HPsiInv <- RootiHPsi %*% t(RootiHPsi)
    #HPsiInv <- chol2inv(CholHPsi)
  }
  
  
  ###Set initial values of Delta
  if ("Delta" %in% UserStarters) {Delta <- matrix(starting$Delta, nrow = K, ncol = 1)} else {Delta <- matrix(1, nrow = K, ncol = 1)} #if ((!"Delta" %in% UserStarters)) 

  ###Set initial values of Beta
  if ("Beta" %in% UserStarters) {Beta <- matrix(starting$Beta, nrow = P, ncol = 1)} else {Beta <- matrix(0, nrow = P, ncol = 1)} #if ((!"Beta" %in% UserStarters)) 
    
  ###Set initial value of Upsilon
  if ("Upsilon" %in% UserStarters) {Upsilon <- matrix(starting$Upsilon, nrow = K, ncol = K)} else {Upsilon <- diag(K)} #if (!("Upsilon" %in% UserStarters)) 
  RootiUpsilon = solve(chol(Upsilon))
  UpsilonInv = RootiUpsilon %*% t(RootiUpsilon)

  ###Create atom variances (Tau2)
  if (GS == 1) {Tau <- matrix(cumprod(Delta), nrow = K, ncol = 1)} else {Tau <- Delta} #if (GS == 0) 
  
  ###Create label parameters
  Xi <- matrix(0, nrow = M * O, ncol = K)
  Lambda <- matrix(0, nrow = M * O, ncol = K)

  ###Factors
  BigPhi <- matrix(0, nrow = K, ncol = Nu)
  Eta <- matrix(as.numeric(BigPhi), ncol = 1)
  
  ###Slice sampling latent parameter
  U <- matrix(runif(M * O * K, 0, 1), nrow = M * O, ncol = K)
 
  ###Spatial covariance objects
  CholKappa <- chol(Kappa)
  RootiKappa <- solve(CholKappa)
  KappaInv <- RootiKappa %*% t(RootiKappa)
  #KappaInv <- chol2inv(CholKappa)

  ###Create Poly-Gamma update objects
  Cov <- array(1, dim = c(M, O, Nu))
  NotCount <- which(FamilyInd != 3)
  count <- 1
  for (f in 1:O) {
    if (FamilyInd[f] %in% 0:2) {
      Cov[ , f, ] <- matrix(Sigma2[, count], nrow = M, ncol = Nu)
      count <- count + 1
    }
    if (FamilyInd[f] == 3) {
      omega <- matrix(pgdraw(rep(1, M*Nu), rep(0, M*Nu)), nrow = M, ncol = Nu)
      Cov[ , f, ] <- 1 / omega
    }
  }

  ###Save parameter objects
  Para <- list()
  Para$Sigma2 <- Sigma2
  Para$Kappa <- Kappa
  Para$Delta <- Delta
  Para$Beta <- Beta
  Para$XBeta <- X %*% Beta
  Para$Psi <- Psi
  Para$Upsilon <- Upsilon
  Para$UpsilonInv <- UpsilonInv
  Para$RootiUpsilon <- RootiUpsilon
  Para$Xi <- Xi
  Para$Lambda <- Lambda
  Para$Tau <- Tau
  Para$BigPhi <- BigPhi
  Para$Eta <- Eta
  Para$HPsi <- HPsi
  Para$RootiHPsi <- RootiHPsi
  Para$HPsiInv <- HPsiInv
  Para$Mean <- kronecker(EyeNu, Lambda) %*% Eta + X %*% Beta
  Para$U <- U
  Para$LjVec <- LjVec
  Para$RootiKappa <- RootiKappa
  Para$KappaInv <- KappaInv
  Para$Cov <- Cov
  return(Para)
}

VAR1CreateParaVaryingLjs <- function(starting, DatObj, HyPara, LjVec) {
  
  ###Set data objects
  K <- DatObj$K
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  FamilyInd <- DatObj$FamilyInd
  X <- DatObj$X
  GS <- DatObj$GS
  
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  ###Set initial values of Sigma2
  if (any(FamilyInd != 3)) {
    if ("Sigma2" %in% UserStarters) {Sigma2 <- matrix(starting$Sigma2, ncol = (O - C), nrow = M)} else {Sigma2 <- matrix(1, ncol = (O - C), nrow = M)}
  } else {Sigma2 <- matrix(1, nrow = 1, ncol = 1)}
  
  ###Set initial values of Kappa
  if ("Kappa" %in% UserStarters) {Kappa <- matrix(starting$Kappa, nrow = O, ncol = O)} else {Kappa <- diag(O)}
  
  ###Set initial values of A
  if ("A" %in% UserStarters) {A <- starting$A} else {A <- 0.8 * diag(K)}
  
  ###Set initial values of Delta
  if ("Delta" %in% UserStarters) {Delta <- matrix(starting$Delta, nrow = K, ncol = 1)} else {Delta <- matrix(1, nrow = K, ncol = 1)}
  
  ###Set initial values of Beta
  if ("Beta" %in% UserStarters) {Beta <- matrix(starting$Beta, nrow = P, ncol = 1)} else {Beta <- matrix(0, nrow = P, ncol = 1)}
  
  ###Set initial value of Upsilon
  if ("Upsilon" %in% UserStarters) {Upsilon <- matrix(starting$Upsilon, nrow = K, ncol = K)} else {Upsilon <- diag(K)}
  RootiUpsilon = solve(chol(Upsilon))
  UpsilonInv = RootiUpsilon %*% t(RootiUpsilon)
  
  ###Create atom variances (Tau2)
  if (GS == 1) {Tau <- matrix(cumprod(Delta), nrow = K, ncol = 1)} else {Tau <- Delta}
  
  ###Create label parameters
  Xi <- matrix(0, nrow = M * O, ncol = K)
  Lambda <- matrix(0, nrow = M * O, ncol = K)
  
  ###Factors
  BigPhi <- matrix(0, nrow = K, ncol = Nu)
  Eta <- matrix(as.numeric(BigPhi), ncol = 1)
  
  ###Slice sampling latent parameter
  U <- matrix(runif(M * O * K, 0, 1), nrow = M * O, ncol = K)
  
  ###Spatial covariance objects
  CholKappa <- chol(Kappa)
  RootiKappa <- solve(CholKappa)
  KappaInv <- RootiKappa %*% t(RootiKappa)
  #KappaInv <- chol2inv(CholKappa)
  
  ###Create Poly-Gamma update objects
  Cov <- array(1, dim = c(M, O, Nu))
  NotCount <- which(FamilyInd != 3)
  count <- 1
  for (f in 1:O) {
    if (FamilyInd[f] %in% 0:2) {
      Cov[ , f, ] <- matrix(Sigma2[, count], nrow = M, ncol = Nu)
      count <- count + 1
    }
    if (FamilyInd[f] == 3) {
      omega <- matrix(pgdraw(rep(1, M*Nu), rep(0, M*Nu)), nrow = M, ncol = Nu)
      Cov[ , f, ] <- 1 / omega
    }
  }
  
  ###Save parameter objects
  Para <- list()
  Para$Sigma2 <- Sigma2
  Para$Kappa <- Kappa
  Para$Delta <- Delta
  Para$Beta <- Beta
  Para$XBeta <- X %*% Beta
  Para$A <- A
  Para$Upsilon <- Upsilon
  Para$UpsilonInv <- UpsilonInv
  Para$Xi <- Xi
  Para$Lambda <- Lambda
  Para$Tau <- Tau
  Para$BigPhi <- BigPhi
  Para$Eta <- Eta
  Para$Mean <- kronecker(EyeNu, Lambda) %*% Eta + X %*% Beta
  Para$U <- U
  Para$LjVec <- LjVec
  Para$RootiKappa <- RootiKappa
  Para$KappaInv <- KappaInv
  Para$Cov <- Cov
  return(Para)
}



## called when !(IS = 1, SpCorInd = 0, SA = 1) for the fixedL main function 
CreateSpatPara1FixedL <- function(starting, DatObj, HyPara){
  ###Set data objects
  M <- DatObj$M
  IS <- DatObj$IS
  SpCorInd <- DatObj$SpCorInd
  SpDist <- DatObj$SpDist
  ###Set hyperparameter objects
  ARho <- HyPara$ARho
  BRho <- HyPara$BRho
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  SpatPara <- list()
  SpatPara$DwVec <- rep(1, M)
  
  ###Set initial values of Rho and Spatial correlation objects
  if (IS == 0){
    EyeM <- diag(M)
    SpCov <- EyeM
    RootiSpCov <- EyeM
    SpCovInv <- EyeM
    Rho <- 0
  } else {# when IS = 1
    if (SpCorInd == 1) { #discrete; spatApprox must be FALSE
      if ("Rho" %in% UserStarters) {Rho <- starting$Rho} else {Rho <- 0.9}
      Dw <- diag(apply(SpDist, 1, sum))
      SpatPara$DwVec <- diag(Dw)
      SpCovInv <- Dw - Rho * SpDist
      RootiSpCov <- t(chol(SpCovInv)) #transpose of the cholesky factor of SpCovInv instead; not the real rooti but works here too
    } else { #SpCorInd = 0 continuous
      if ("Rho" %in% UserStarters){
        Rho <- starting$Rho
        if ((Rho <= ARho) | (Rho >= BRho)) stop('starting: "Rho" must be in (ARho, BRho)')
      } else {
        Rho <- mean(c(ARho, BRho))
      }
      SpCov <- exp(-Rho * SpDist)
      ## since SA must = 0 if this function CreateSpatPara1 is called and IS = 1, SpCorInd = 0
      RootiSpCov <- solve(chol(SpCov))
      SpCovInv <- RootiSpCov %*% t(RootiSpCov)
    }
  }
  
  ###Save parameter objects
  SpatPara$Rho <- Rho
  SpatPara$RootiSpCov <- RootiSpCov
  SpatPara$SpCovInv <- SpCovInv
  return(SpatPara)
}

## called when !(IS = 1, SpCorInd = 0, SA = 1) for the VaryingLjs main function 
CreateSpatPara1VaryingLjs <- function(starting, DatObj, HyPara){
  ###Set data objects
  M <- DatObj$M
  IS <- DatObj$IS
  SpCorInd <- DatObj$SpCorInd
  SpDist <- DatObj$SpDist
  ###Set hyperparameter objects
  ARho <- HyPara$ARho
  BRho <- HyPara$BRho
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  SpatPara <- list()
  SpatPara$DwVec <- rep(1, M)
  
  ###Set initial values of Rho and Spatial correlation objects
  if (IS == 0){
    EyeM <- diag(M)
    SpCov <- EyeM
    RootiSpCov <- EyeM
    SpCovInv <- EyeM
    Rho <- 0
  } else {# when IS = 1
    if (SpCorInd == 1) { #discrete; spatApprox must be FALSE
      if ("Rho" %in% UserStarters) {Rho <- starting$Rho} else {Rho <- 0.9}
      Dw <- diag(apply(SpDist, 1, sum))
      SpatPara$DwVec <- diag(Dw)
      SpCovInv <- Dw - Rho * SpDist
      SpCov <- solve(SpCovInv)
      RootiSpCov <- t(chol(SpCovInv)) #transpose of the cholesky factor of SpCovInv instead; not the real rooti but works here too
    } else { #SpCorInd = 0 continuous
      if ("Rho" %in% UserStarters){
        Rho <- starting$Rho
        if ((Rho <= ARho) | (Rho >= BRho)) stop('starting: "Rho" must be in (ARho, BRho)')
      } else {
        Rho <- mean(c(ARho, BRho))
      }
      SpCov <- exp(-Rho * SpDist)
      ## since SA must = 0 if this function CreateSpatPara1 is called and IS = 1, SpCorInd = 0
      RootiSpCov <- solve(chol(SpCov))
      SpCovInv <- RootiSpCov %*% t(RootiSpCov)
    }
  }
  
  ###Save parameter objects
  SpatPara$Rho <- Rho
  SpatPara$SpCov <- SpCov
  SpatPara$RootiSpCov <- RootiSpCov
  SpatPara$SpCovInv <- SpCovInv
  return(SpatPara)
}


## called when IS = 1, SpCorInd = 0, SA = 1, and (CL = 0 or (CL = 1 and alphaMethod = "block")) for the fixedL main function 
CreateSpatPara2FixedL <- function(starting, DatObj, HyPara){
  ###Set data objects
  M <- DatObj$M
  SpDist <- DatObj$SpDist
  h <- DatObj$h
  ###Set hyperparameter objects
  ARho <- HyPara$ARho
  BRho <- HyPara$BRho
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  ###Set initial values of Rho and Spatial correlation objects
  if ("Rho" %in% UserStarters){
    Rho <- starting$Rho
    if ((Rho <= ARho) | (Rho >= BRho)) stop('starting: "Rho" must be in (ARho, BRho)')
  } else {
    Rho <- mean(c(ARho, BRho))
  }
  SpCov <- exp(-Rho * SpDist)
  nnInd <- matrix(0, M, h)
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nnInd[i, 1:nnMaxNum] <- sort(order(SpDist[i, 1:i-1])[1:nnMaxNum])
  }
  Bs <- matrix(0, M, h)
  Fs <- vector()
  Fs[1] = 1
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nni <- nnInd[i, 1:nnMaxNum]
    Bs[i, 1:nnMaxNum] <- SpCov[i, nni]%*%solve(SpCov[nni, nni])
    Fs[i] <- 1 - Bs[i, 1:nnMaxNum]%*%SpCov[nni, i]
  }
  detSpCovApprox <- prod(Fs)
  BsStar <- matrix(0, M, M)
  diag(BsStar) <- rep(1, M)
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nni <- nnInd[i, 1:nnMaxNum]
    BsStar[i, nni] <- -Bs[i, 1:nnMaxNum]
  }
  SpCovInvApprox <- t(BsStar)%*%diag(Fs^{-1})%*%BsStar
  
  ###Save parameter objects
  SpatPara <- list()
  SpatPara$Rho <- Rho
  SpatPara$nnInd <- nnInd - 1
  SpatPara$detSpCov <- detSpCovApprox
  SpatPara$SpCovInv <- SpCovInvApprox
  return(SpatPara)
} 

## called when IS = 1, SpCorInd = 0, SA = 1, and alphaSequen = FALSE for the VaryingLjs main function
CreateSpatPara2VaryingLjs <- function(starting, DatObj, HyPara){
  ###Set data objects
  M <- DatObj$M
  SpDist <- DatObj$SpDist
  h <- DatObj$h
  ###Set hyperparameter objects
  ARho <- HyPara$ARho
  BRho <- HyPara$BRho
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  ###Set initial values of Rho and Spatial correlation objects
  if ("Rho" %in% UserStarters){
    Rho <- starting$Rho
    if ((Rho <= ARho) | (Rho >= BRho)) stop('starting: "Rho" must be in (ARho, BRho)')
  } else {
    Rho <- mean(c(ARho, BRho))
  }
  SpCov <- exp(-Rho * SpDist)
  nnInd <- matrix(0, M, h)
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nnInd[i, 1:nnMaxNum] <- sort(order(SpDist[i, 1:i-1])[1:nnMaxNum])
  }
  Bs <- matrix(0, M, h)
  Fs <- vector()
  Fs[1] = 1
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nni <- nnInd[i, 1:nnMaxNum]
    Bs[i, 1:nnMaxNum] <- SpCov[i, nni]%*%solve(SpCov[nni, nni])
    Fs[i] <- 1 - Bs[i, 1:nnMaxNum]%*%SpCov[nni, i]
  }
  detSpCovApprox <- prod(Fs)
  BsStar <- matrix(0, M, M)
  diag(BsStar) <- rep(1, M)
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nni <- nnInd[i, 1:nnMaxNum]
    BsStar[i, nni] <- -Bs[i, 1:nnMaxNum]
  }
  SpCovInvApprox <- t(BsStar)%*%diag(Fs^{-1})%*%BsStar
  
  ###Save parameter objects
  SpatPara <- list()
  SpatPara$Rho <- Rho
  SpatPara$SpCov <- SpCov
  SpatPara$nnInd <- nnInd - 1
  SpatPara$detSpCov <- detSpCovApprox
  SpatPara$SpCovInv <- SpCovInvApprox
  return(SpatPara)
} 


## called when IS = 1, SpCorInd = 0, SA = 1, CL = 1, and alphaMethod = "sequential" for the fixedL main function and
# when IS = 1, SpCorInd = 0, SA = 1, and alphaSequen = TRUE for the VaryingLjs main function 
CreateSpatPara3 <- function(starting, DatObj, HyPara){
  ###Set data objects
  M <- DatObj$M
  SpDist <- DatObj$SpDist
  h <- DatObj$h
  ###Set hyperparameter objects
  ARho <- HyPara$ARho
  BRho <- HyPara$BRho
  ###Which parameters are user defined?
  UserStarters <- names(starting)
  
  ###Set initial values of Rho and Spatial correlation objects
  if ("Rho" %in% UserStarters){
    Rho <- starting$Rho
    if ((Rho <= ARho) | (Rho >= BRho)) stop('starting: "Rho" must be in (ARho, BRho)')
  } else {
    Rho <- mean(c(ARho, BRho))
  }
  SpCov <- exp(-Rho * SpDist)
  nnInd <- matrix(0, M, h)
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nnInd[i, 1:nnMaxNum] <- sort(order(SpDist[i, 1:i-1])[1:nnMaxNum])
  }
  Bs <- matrix(0, M, h)
  Fs <- vector()
  Fs[1] = 1
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nni <- nnInd[i, 1:nnMaxNum]
    Bs[i, 1:nnMaxNum] <- SpCov[i, nni]%*%solve(SpCov[nni, nni])
    Fs[i] <- 1 - Bs[i, 1:nnMaxNum]%*%SpCov[nni, i]
  }
  detSpCovApprox <- prod(Fs)
  BsStar <- matrix(0, M, M)
  diag(BsStar) <- rep(1, M)
  for(i in 2:M){
    nnMaxNum <- min(h, i-1)
    nni <- nnInd[i, 1:nnMaxNum]
    BsStar[i, nni] <- -Bs[i, 1:nnMaxNum]
  }
  SpCovInvApprox <- t(BsStar)%*%diag(Fs^{-1})%*%BsStar
  
  ###Save parameter objects
  SpatPara <- list()
  SpatPara$Rho <- Rho
  SpatPara$Bs <- Bs
  SpatPara$Fs <- Fs
  SpatPara$nnInd <- nnInd - 1
  SpatPara$detSpCov <- detSpCovApprox
  SpatPara$SpCovInv <- SpCovInvApprox
  return(SpatPara)
}



###Function that creates the data augmentation (i.e. Tobit) booleans------------------------------------------------------------
CreateDatAug <- function(DatObj) {

  ###Set data object
  YObserved <- DatObj$YObserved
  FamilyInd <- DatObj$FamilyInd
  YStarWide <- DatObj$YStarWide
  M <- DatObj$M
  O <- DatObj$O
  Nu <- DatObj$Nu
  N <- DatObj$N

  ###Initialize Data Augmentation Object with Normality
  DatAug <- list()
  DatAug$NBelow <- 0
  DatAug$NAbove <- 0
  DatAug$WhichBelow <- 0
  DatAug$WhichAbove <- 0
  
  ###Lower truncation censored from below at 0
  # the recorded observed Y is 0 (0 is the truncated value) is equivalent to 
  # the actual pre-truncated observed Y value < 0 
  if (any(FamilyInd == 2) | any(FamilyInd == 1)) {#Tobit | Probit 
    TobitBoolean <- matrix(FALSE, ncol = 1, nrow = N)
    for (o in 1:O) {
      if (FamilyInd[o] == 2) {
        TobitBooleanCube <- array(FALSE, dim = c(M, O, Nu))
        TobitBooleanCube[ , o, ] <- TRUE
        TobitBoolean <- TobitBoolean | matrix(TobitBooleanCube & YObserved == 0, ncol = 1)
      }
    }
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    DatAug$NBelow <- NBelow
    DatAug$WhichBelow <- WhichBelow - 1
  }
  
  ###Upper truncation Probit
  # the recorded observed probit Y is 1 is equivalent to 
  # the actual pre-truncated observed Y value > 0 
  if (any(FamilyInd == 1)) {#Probit
    ProbitBoolean <- matrix(FALSE, ncol = 1, nrow = N)
    for (o in 1:O) {
      if (FamilyInd[o] == 1) {
        ProbitBooleanCube <- array(FALSE, dim = c(M, O, Nu))
        ProbitBooleanCube[ , o, ] <- TRUE
        ProbitBoolean <- ProbitBoolean | matrix(ProbitBooleanCube & YObserved == 1, ncol = 1)
      }
    }
    WhichAbove <- which(ProbitBoolean)
    NAbove <- length(WhichAbove)
    DatAug$WhichAbove <- WhichAbove - 1
    DatAug$NAbove <- NAbove
  }
  return(DatAug)

}



###Function that creates inputs for MCMC sampler--------------------------------------------------------------------------------
CreateMcmc <- function(mcmc, DatObj) {

  ###Which parameters are user defined?
  UserMCMC <- names(mcmc)

  ###Set MCMC objects
  if ("NBurn" %in% UserMCMC) NBurn <- mcmc$NBurn else NBurn <- 10000
  if ("NSims" %in% UserMCMC) NSims <- mcmc$NSims else NSims <- 10000
  if ("NThin" %in% UserMCMC) NThin <- mcmc$NThin else NThin <- 1
  if ("NPilot" %in% UserMCMC) NPilot <- mcmc$NPilot else NPilot <- 20

  ###One last check of MCMC user inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!(is.wholenumber(NSims / NThin))) stop('mcmc: "NThin" must be a factor of "NSims"')
  if (!(is.wholenumber(NBurn / NPilot))) stop('mcmc: "NPilot" must be a factor of "NBurn"')

  ###Create MCMC objects
  NTotal <- NBurn + NSims
  WhichKeep <- NBurn + (1:(NSims / NThin)) * NThin
  NKeep <- length(WhichKeep)

  ###Pilot adaptation objects
  WhichPilotAdapt <- (1:NPilot) * NBurn / NPilot
  PilotAdaptDenominator <- WhichPilotAdapt[1]

  ###Burn-in progres bar
  BarLength <- 100 #Burn-in bar length (arbitrary)
  BurnInProgress <- seq(1 / BarLength, 1, 1 / BarLength)
  WhichBurnInProgress <- sapply(BurnInProgress, function(x) tail(which(1 : NBurn <= x * NBurn), 1))

  ###Progress output objects
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichSamplerProgress <- sapply(SamplerProgress, function(x) tail(which(1:NSims <= x * NSims), 1)) + NBurn
  WhichBurnInProgressInt <- sapply(SamplerProgress, function(x) tail(which(1:NBurn <= x * NBurn), 1))

  ###Save objects
  MCMC <- list()
  MCMC$NBurn <- NBurn
  MCMC$NSims <- NSims
  MCMC$NThin <- NThin
  MCMC$NPilot <- NPilot
  MCMC$NTotal <- NTotal
  MCMC$WhichKeep <- WhichKeep
  MCMC$NKeep <- NKeep
  MCMC$WhichPilotAdapt <- WhichPilotAdapt
  MCMC$PilotAdaptDenominator <- PilotAdaptDenominator
  MCMC$BurnInProgress <- BurnInProgress
  MCMC$WhichBurnInProgress <- WhichBurnInProgress
  MCMC$WhichBurnInProgressInt <- WhichBurnInProgressInt
  MCMC$BarLength <- BarLength
  MCMC$WhichSamplerProgress <- WhichSamplerProgress
  return(MCMC)

}



###Functions that create memory spaces for storing posterior estimates for parameters with fixed dimensions -------------------------------------------------
CreateStorageFixedL <- function(DatObj, McmcObj) {
  
  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  L <- DatObj$L
  
  ###Set MCMC objects
  NKeep <- McmcObj$NKeep
  
  ###Create storage object 
  # Lambda: M * O * K
  # Eta: K * Nu
  # Sigma: M * (O - C)
  # Kappa: O * (O + 1) / 2
  # Rho: 1
  # Upsilon: K * (K + 1) / 2
  # Psi: 1
  # Beta: P
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + 1 + 1 + P
  if (DatObj$CL == 1){
    # Xi: M * O * K
    # Delta: K
    # Alpha: M * O * K * (L - 1)
    # Theta: K * L
    if (DatObj$spatPred == 1) {
      numPara = numPara + M * O * K + K + K * L
      if (DatObj$alphaWsFiles == 0) numPara = numPara + M * O * K * (L - 1) 
    }
    # Weights: M * O * K * L
    if ((DatObj$storeW == 1) & (DatObj$alphaWsFiles == 0)) numPara = numPara + M * O * K * L
  }
  
  Out <- matrix(nrow = numPara, ncol = NKeep)
  return(Out)
  
}

VAR1CreateStorageFixedL <- function(DatObj, McmcObj) {
  
  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  L <- DatObj$L
  
  ###Set MCMC objects
  NKeep <- McmcObj$NKeep
  
  ###Create storage object 
  # Lambda: M * O * K
  # Eta: K * Nu
  # Sigma: M * (O - C)
  # Kappa: O * (O + 1) / 2
  # Rho: 1
  # Upsilon: K * (K + 1) / 2
  # A: K * K
  # Beta: P
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + 1 + K^2 + P
  if (DatObj$CL == 1){
    # Xi: M * O * K
    # Delta: K
    # Alpha: M * O * K * (L - 1)
    # Theta: K * L
    if (DatObj$spatPred == 1) {
      numPara = numPara + M * O * K + K + K * L
      if (DatObj$alphaWsFiles == 0) numPara = numPara + M * O * K * (L - 1) 
    }
    # Weights: M * O * K * L
    if ((DatObj$storeW == 1) & (DatObj$alphaWsFiles == 0)) numPara = numPara + M * O * K * L
  }
  
  Out <- matrix(nrow = numPara, ncol = NKeep)
  return(Out)
  
}

CreateStorageVaryingLjs <- function(DatObj, McmcObj) {
  
  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  
  ###Set MCMC objects
  NKeep <- McmcObj$NKeep
  
  ###Create storage object 
  # Lambda: M * O * K
  # Eta: K * Nu
  # Sigma: M * (O - C)
  # Kappa: O * (O + 1) / 2
  # Rho: 1
  # Upsilon: K * (K + 1) / 2
  # Psi: 1
  # Beta: P
  # Xi: M * O * K
  # Delta: K
  # LjVec: K
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + 1 + 1 + P + M * O * K + K + K
  Out <- matrix(nrow = numPara, ncol = NKeep)
  return(Out)
  
}

VAR1CreateStorageVaryingLjs <- function(DatObj, McmcObj) {
  
  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  P <- DatObj$P
  C <- DatObj$C
  
  ###Set MCMC objects
  NKeep <- McmcObj$NKeep
  
  ###Create storage object 
  # Lambda: M * O * K
  # Eta: K * Nu
  # Sigma: M * (O - C)
  # Kappa: O * (O + 1) / 2
  # Rho: 1
  # Upsilon: K * (K + 1) / 2
  # A: K * K
  # Beta: P
  # Xi: M * O * K
  # Delta: K
  # LjVec: K
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + 1 + K^2 + P + M * O * K + K + K
  Out <- matrix(nrow = numPara, ncol = NKeep)
  return(Out)
  
}