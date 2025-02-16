###Functions for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamplesFixedL <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  C <- DatObj$C
  P <- DatObj$P
  L <- DatObj$L
  GS <- DatObj$GS
  CL <- DatObj$CL
  spatPred <- DatObj$spatPred
  storeWeights <- DatObj$storeW
  alphasWeightsToFiles <- DatObj$alphaWsFiles
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(O * M * K)]
  Eta <- RawSamples[, (O * M * K + 1):(O * M * K + K * Nu), drop = FALSE]
  if (C != O) { #if (C == O) 
    Sigma2 <- RawSamples[, (O * M * K + K * Nu + 1):(O * M * K + K * Nu + M * (O - C))]
    Sigma2Ind <- expand.grid(which(DatObj$FamilyInd != 3), 1:M)
    colnames(Sigma2) <- paste0("Sigma2_", Sigma2Ind[, 1], "_", Sigma2Ind[, 2])
  } else {
    Sigma2 <- NULL
  }
  Kappa <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2), drop = FALSE]
  Upsilon <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2), drop = FALSE]
  Psi <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1), drop = FALSE]
  Rho <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1 + 1), drop = FALSE]
  if (P > 0) { #if (P == 0)
    Beta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1 + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1 + 1 + P), drop = FALSE]
    colnames(Beta) <- paste0("Beta", 1:P)
  } else {
    Beta <- NULL
  }
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + 1 + 1 + P
  if (CL == 1) {
    if (spatPred == 1) {
      Xi <- RawSamples[, (numPara + 1):(numPara + O * M * K)]
      XiInd <- expand.grid(1:K, 1:M, 1:O)
      colnames(Xi) <- paste0("Xi_", XiInd[, 3], "_", XiInd[, 2], "_", XiInd[, 1])
      Delta <- RawSamples[, (numPara + O * M * K + 1):(numPara + O * M * K + K), drop = FALSE]
      colnames(Delta) <- paste0("Delta", 1:K)
      if (GS == 1) {Tau <- matrix(t(apply(Delta, 1, cumprod)), ncol = K)} else {Tau <- Delta}
      colnames(Tau) <- paste0("Tau", 1:K)
      if (alphasWeightsToFiles == 0){
        Alpha <- RawSamples[, (numPara + O * M * K + K + 1):(numPara + O * M * K + K + M * O * K * (L - 1)), drop = FALSE]
        Theta <- RawSamples[, (numPara + O * M * K + K + M * O * K * (L - 1) + 1):(numPara + O * M * K + K + M * O * K * (L - 1) + K * L), drop = FALSE]
        AlphaInd <- expand.grid(1:O, 1:M, 1:(L-1), 1:K)
        colnames(Alpha) <- paste("Alpha", AlphaInd[, 4], AlphaInd[, 3], AlphaInd[, 2], AlphaInd[, 1], sep = "_")
        ThetaInd <- expand.grid(1:L, 1:K)
        colnames(Theta) <- paste("Theta", ThetaInd[, 2], ThetaInd[, 1], sep = "_")
        numPara = numPara + M * O * K + K + M * O * K * (L - 1) + K * L
      } else {
        Alpha <- NULL
        Theta <- RawSamples[, (numPara + O * M * K + K + 1):(numPara + O * M * K + K + K * L), drop = FALSE]
        ThetaInd <- expand.grid(1:L, 1:K)
        colnames(Theta) <- paste("Theta", ThetaInd[, 2], ThetaInd[, 1], sep = "_")
        numPara = numPara + M * O * K + K + K * L
      }
    } else {
      Tau <- Delta <- Xi <- Alpha <- Theta <- NULL
    }
    if ((storeWeights == 1) & (alphasWeightsToFiles == 0)){
      Weights <- RawSamples[, (numPara + 1):(numPara + M * O * K * L), drop = FALSE]
      WeightsInd <- expand.grid(1:M, 1:O, 1:L, 1:K)
      colnames(Weights) <- paste("Weights", WeightsInd[, 4], WeightsInd[, 3], WeightsInd[, 2], WeightsInd[, 1], sep = "_")
    } else {
      Weights <- NULL
    }
  } else{
    Tau <- Delta <- Xi <- Alpha <- Theta <- Weights <- NULL
  }
  
  LambdaInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Lambda) <- paste0("Lambda_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  EtaInd <- expand.grid(1:K, 1:Nu)
  colnames(Eta) <- paste0("Eta", EtaInd[, 2], "_", EtaInd[, 1])
  KappaInd <- which(lower.tri(apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (O == 1) {
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]][1]
  } else{ #if (O > 1) 
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]]
  }
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (K == 1) {
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]][1]
  } else { #if (K > 1) 
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  }
  colnames(Psi) <- "Psi"
  colnames(Rho) <- "Rho"
  
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa = Kappa, Delta = Delta, Tau = Tau, Upsilon = Upsilon, Psi = Psi, Xi = Xi, Rho = Rho, 
              Beta = Beta, Theta = Theta, Alpha = Alpha, Weights = Weights)
  return(Out)
}

VAR1FormatSamplesFixedL <- function(DatObj, RawSamples) {
  
  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  C <- DatObj$C
  P <- DatObj$P
  L <- DatObj$L
  GS <- DatObj$GS
  CL <- DatObj$CL
  spatPred <- DatObj$spatPred
  storeWeights <- DatObj$storeW
  alphasWeightsToFiles <- DatObj$alphaWsFiles
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(O * M * K)]
  Eta <- RawSamples[, (O * M * K + 1):(O * M * K + K * Nu), drop = FALSE]
  if (C != O) { #if (C == O) 
    Sigma2 <- RawSamples[, (O * M * K + K * Nu + 1):(O * M * K + K * Nu + M * (O - C))]
    Sigma2Ind <- expand.grid(which(DatObj$FamilyInd != 3), 1:M)
    colnames(Sigma2) <- paste0("Sigma2_", Sigma2Ind[, 1], "_", Sigma2Ind[, 2])
  } else {
    Sigma2 <- NULL
  }
  Kappa <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2), drop = FALSE]
  Upsilon <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2), drop = FALSE]
  A <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2), drop = FALSE]
  Rho <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2 + 1), drop = FALSE]
  if (P > 0) { #if (P == 0)
    Beta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2 + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2 + 1 + P), drop = FALSE]
    colnames(Beta) <- paste0("Beta", 1:P)
  } else {
    Beta <- NULL
  }
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + K^2 + 1 + P
  if (CL == 1) {
    if (spatPred == 1) {
      Xi <- RawSamples[, (numPara + 1):(numPara + O * M * K)]
      XiInd <- expand.grid(1:K, 1:M, 1:O)
      colnames(Xi) <- paste0("Xi_", XiInd[, 3], "_", XiInd[, 2], "_", XiInd[, 1])
      Delta <- RawSamples[, (numPara + O * M * K + 1):(numPara + O * M * K + K), drop = FALSE]
      colnames(Delta) <- paste0("Delta", 1:K)
      if (GS == 1) {Tau <- matrix(t(apply(Delta, 1, cumprod)), ncol = K)} else {Tau <- Delta}
      colnames(Tau) <- paste0("Tau", 1:K)
      if (alphasWeightsToFiles == 0){
        Alpha <- RawSamples[, (numPara + O * M * K + K + 1):(numPara + O * M * K + K + M * O * K * (L - 1)), drop = FALSE]
        Theta <- RawSamples[, (numPara + O * M * K + K + M * O * K * (L - 1) + 1):(numPara + O * M * K + K + M * O * K * (L - 1) + K * L), drop = FALSE]
        AlphaInd <- expand.grid(1:O, 1:M, 1:(L-1), 1:K)
        colnames(Alpha) <- paste("Alpha", AlphaInd[, 4], AlphaInd[, 3], AlphaInd[, 2], AlphaInd[, 1], sep = "_")
        ThetaInd <- expand.grid(1:L, 1:K)
        colnames(Theta) <- paste("Theta", ThetaInd[, 2], ThetaInd[, 1], sep = "_")
        numPara = numPara + M * O * K + K + M * O * K * (L - 1) + K * L
      } else {
        Alpha <- NULL
        Theta <- RawSamples[, (numPara + O * M * K + K + 1):(numPara + O * M * K + K + K * L), drop = FALSE]
        ThetaInd <- expand.grid(1:L, 1:K)
        colnames(Theta) <- paste("Theta", ThetaInd[, 2], ThetaInd[, 1], sep = "_")
        numPara = numPara + M * O * K + K + K * L
      }
    } else {
      Tau <- Delta <- Xi <- Alpha <- Theta <- NULL
    }
    if ((storeWeights == 1) & (alphasWeightsToFiles == 0)){
      Weights <- RawSamples[, (numPara + 1):(numPara + M * O * K * L), drop = FALSE]
      WeightsInd <- expand.grid(1:M, 1:O, 1:L, 1:K)
      colnames(Weights) <- paste("Weights", WeightsInd[, 4], WeightsInd[, 3], WeightsInd[, 2], WeightsInd[, 1], sep = "_")
    } else {
      Weights <- NULL
    }
  } else{
    Tau <- Delta <- Xi <- Alpha <- Theta <- Weights <- NULL
  }
  
  LambdaInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Lambda) <- paste0("Lambda_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  EtaInd <- expand.grid(1:K, 1:Nu)
  colnames(Eta) <- paste0("Eta", EtaInd[, 2], "_", EtaInd[, 1])
  KappaInd <- which(lower.tri(apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (O == 1) {
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]][1]
  } else{ #if (O > 1) 
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]]
  }
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (K == 1) {
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]][1]
  } else { #if (K > 1) 
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  }
  AInd <- expand.grid(1:K, 1:K)
  colnames(A) <- paste0("A", "_row", AInd[, 2], "_col",AInd[, 1])
  colnames(Rho) <- "Rho"
  
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa = Kappa, Delta = Delta, Tau = Tau, Upsilon = Upsilon, A = A, Xi = Xi, Rho = Rho, 
              Beta = Beta, Theta = Theta, Alpha = Alpha, Weights = Weights)
  return(Out)
}


FormatSamplesVaryingLjs <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  C <- DatObj$C
  P <- DatObj$P
  GS <- DatObj$GS
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(O * M * K)]
  Eta <- RawSamples[, (O * M * K + 1):(O * M * K + K * Nu), drop = FALSE]
  if (C != O) {
    Sigma2 <- RawSamples[, (O * M * K + K * Nu + 1):(O * M * K + K * Nu + M * (O - C))]
    Sigma2Ind <- expand.grid(which(DatObj$FamilyInd != 3), 1:M)
    colnames(Sigma2) <- paste0("Sigma2_", Sigma2Ind[, 1], "_", Sigma2Ind[, 2])
  } else { # if (C == O) 
    Sigma2 <- NULL
  }
  Kappa <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2), drop = FALSE]
  Upsilon <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2), drop = FALSE]
  Psi <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1), drop = FALSE]
  Rho <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1 + 1), drop = FALSE]
  if (P > 0) {
    Beta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1 + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1 + 1 + P), drop = FALSE]
    colnames(Beta) <- paste0("Beta", 1:P)
  } else { # if (P == 0) 
    Beta <- NULL
  }
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + 1 + 1 + P
  Xi <- RawSamples[, (numPara + 1):(numPara + O * M * K)]
  XiInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Xi) <- paste0("Xi_", XiInd[, 3], "_", XiInd[, 2], "_", XiInd[, 1])
  Delta <- RawSamples[, (numPara + O * M * K + 1):(numPara + O * M * K + K), drop = FALSE]
  colnames(Delta) <- paste0("Delta", 1:K)
  if (GS == 1) {Tau <- matrix(t(apply(Delta, 1, cumprod)), ncol = K)} else {Tau <- Delta}
  colnames(Tau) <- paste0("Tau", 1:K)
  LjVec <- RawSamples[, (numPara + O * M * K + K + 1):(numPara + O * M * K + K + K), drop = FALSE]
  colnames(LjVec) <- paste0("LjVec", 1:K)
  numPara = numPara + M * O * K + K + K
  
  LambdaInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Lambda) <- paste0("Lambda_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  EtaInd <- expand.grid(1:K, 1:Nu)
  colnames(Eta) <- paste0("Eta", EtaInd[, 2], "_", EtaInd[, 1])
  KappaInd <- which(lower.tri(apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (O == 1) {
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]][1]
  } else { #if (O > 1) 
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]]
  }
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (K == 1) {
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]][1]
  } else {
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  }
  colnames(Psi) <- "Psi"
  colnames(Rho) <- "Rho"
  
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa = Kappa, Delta = Delta, Tau = Tau, Upsilon = Upsilon, Psi = Psi, Xi = Xi, Rho = Rho, Beta = Beta, LjVec = LjVec)
  return(Out)
}

VAR1FormatSamplesVaryingLjs <- function(DatObj, RawSamples) {
  
  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  C <- DatObj$C
  P <- DatObj$P
  GS <- DatObj$GS
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(O * M * K)]
  Eta <- RawSamples[, (O * M * K + 1):(O * M * K + K * Nu), drop = FALSE]
  if (C != O) {
    Sigma2 <- RawSamples[, (O * M * K + K * Nu + 1):(O * M * K + K * Nu + M * (O - C))]
    Sigma2Ind <- expand.grid(which(DatObj$FamilyInd != 3), 1:M)
    colnames(Sigma2) <- paste0("Sigma2_", Sigma2Ind[, 1], "_", Sigma2Ind[, 2])
  } else { # if (C == O) 
    Sigma2 <- NULL
  }
  Kappa <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2), drop = FALSE]
  Upsilon <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2), drop = FALSE]
  A <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2), drop = FALSE]
  Rho <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2 + 1), drop = FALSE]
  if (P > 0) {
    Beta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2 + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + (K * (K + 1)) / 2 + K^2 + 1 + P), drop = FALSE]
    colnames(Beta) <- paste0("Beta", 1:P)
  } else { # if (P == 0) 
    Beta <- NULL
  }
  numPara = M * O * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + ((K + 1) * K) / 2 + K^2 + 1 + P
  Xi <- RawSamples[, (numPara + 1):(numPara + O * M * K)]
  XiInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Xi) <- paste0("Xi_", XiInd[, 3], "_", XiInd[, 2], "_", XiInd[, 1])
  Delta <- RawSamples[, (numPara + O * M * K + 1):(numPara + O * M * K + K), drop = FALSE]
  colnames(Delta) <- paste0("Delta", 1:K)
  if (GS == 1) {Tau <- matrix(t(apply(Delta, 1, cumprod)), ncol = K)} else {Tau <- Delta}
  colnames(Tau) <- paste0("Tau", 1:K)
  LjVec <- RawSamples[, (numPara + O * M * K + K + 1):(numPara + O * M * K + K + K), drop = FALSE]
  colnames(LjVec) <- paste0("LjVec", 1:K)
  numPara = numPara + M * O * K + K + K
  
  LambdaInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Lambda) <- paste0("Lambda_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  EtaInd <- expand.grid(1:K, 1:Nu)
  colnames(Eta) <- paste0("Eta", EtaInd[, 2], "_", EtaInd[, 1])
  KappaInd <- which(lower.tri(apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (O == 1) {
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]][1]
  } else { #if (O > 1) 
    colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]]
  }
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (K == 1) {
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]][1]
  } else {
    colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  }
  AInd <- expand.grid(1:K, 1:K)
  colnames(A) <- paste0("A", "_row", AInd[, 2], "_col",AInd[, 1])
  colnames(Rho) <- "Rho"
  
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa = Kappa, Delta = Delta, Tau = Tau, Upsilon = Upsilon, A = A, Xi = Xi, Rho = Rho, Beta = Beta, LjVec = LjVec)
  return(Out)
}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  IS <- DatObj$IS

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropPsi <- MetropRcpp$MetropPsi
  AcceptancePsi <- MetropRcpp$AcceptancePsi
  if (IS == 0){
    OriginalTuners <- MetrObj$OriginalTuners[1]
    AcceptancePct <- AcceptancePsi / NSims
    MetrSummary <- cbind(AcceptancePct, MetropPsi, OriginalTuners)
    rownames(MetrSummary) <- "Psi"
    colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  } else{ #add rho if IS == 1
    MetropRho <- MetropRcpp$MetropRho
    AcceptanceRho <- MetropRcpp$AcceptanceRho
    OriginalTuners <- MetrObj$OriginalTuners
    AcceptancePct <- c(AcceptancePsi, AcceptanceRho) / NSims
    MetrSummary <- cbind(AcceptancePct, c(MetropPsi, MetropRho), OriginalTuners)
    rownames(MetrSummary) <- c("Psi", "Rho")
    colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  }

  ###Summarize and output
  return(MetrSummary)

}

VAR1SummarizeMetropolis <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {
  
    NSims <- McmcObj$NSims
    MetropRho <- MetropRcpp$MetropRho
    AcceptanceRho <- MetropRcpp$AcceptanceRho
    OriginalTuner <- MetrObj$OriginalTuner
    AcceptancePct <- AcceptanceRho / NSims
    MetrSummary <- cbind(AcceptancePct, MetropRho, OriginalTuner)
    rownames(MetrSummary) <- "Rho"
    colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuner", "OriginalTuner")
    return(MetrSummary)
  
}







###Verify the class of our regression object------------------------------------------------------------------------



#' is.FixedLbfa
#'
#' \code{is.FixedLbfa} is a general test of an object being interpretable as a
#' \code{\link{FixedLbfa}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{FixedLbfa}} class is defined as the regression object that
#'  results from the \code{\link{bfaFixedL}} regression function.
#'
#' @return \code{is.FixedLbfa} returns a logical, depending on whether the object is of class \code{\link{FixedLbfa}}.
#'   
#' @examples 
#' ###Load pre-computed results
#' data(reg.FixedLbfa)
#' 
#' ###Test function
#' is.FixedLbfa(reg.FixedLbfa)
#'
#' @export
is.FixedLbfa <- function(x) {
  identical(attributes(x)$class, "FixedLbfa")
}

#' is.varyLjBFA
#'
#' \code{is.varyLjBFA} is a general test of an object being interpretable as a
#' \code{\link{varyLjBFA}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{varyLjBFA}} class is defined as the regression object that
#'  results from the \code{\link{bfaVaryingLjs}} regression function.
#'
#' @return \code{is.varyLjBFA} returns a logical, depending on whether the object is of class \code{\link{varyLjBFA}}.
#'   
#' @examples 
#' ###Load pre-computed results
#' data(reg.varyLjBFA)
#' 
#' ###Test function
#' is.varyLjBFA(reg.varyLjBFA)
#'
#' @export
is.varyLjBFA <- function(x) {
  identical(attributes(x)$class, "varyLjBFA")
}



#' is.FixedLbfaVAR1
#'
#' \code{is.FixedLbfaVAR1} is a general test of an object being interpretable as a
#' \code{\link{FixedLbfaVAR1}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{FixedLbfaVAR1}} class is defined as the regression object that
#'  results from the \code{\link{VAR1bfaFixedL}} regression function.
#'
#' @return \code{is.FixedLbfaVAR1} returns a logical, depending on whether the object is of class \code{\link{FixedLbfaVAR1}}.
#'   
#' @examples 
#' ###Load pre-computed results
#' data(reg.FixedLbfaVAR1)
#' 
#' ###Test function
#' is.FixedLbfaVAR1(reg.FixedLbfaVAR1)
#'
#' @export
is.FixedLbfaVAR1 <- function(x) {
  identical(attributes(x)$class, "FixedLbfaVAR1")
}

#' is.VAR1varyLjBFA
#'
#' \code{is.VAR1varyLjBFA} is a general test of an object being interpretable as a
#' \code{\link{VAR1varyLjBFA}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{VAR1varyLjBFA}} class is defined as the regression object that
#'  results from the \code{\link{VAR1bfaVaryingLjs}} regression function.
#'
#' @return \code{is.VAR1varyLjBFA} returns a logical, depending on whether the object is of class \code{\link{VAR1varyLjBFA}}.
#'   
#' @examples 
#' ###Load pre-computed results
#' data(reg.VAR1varyLjBFA)
#' 
#' ###Test function
#' is.varyLjBFA(reg.VAR1varyLjBFA)
#'
#' @export
is.VAR1varyLjBFA <- function(x) {
  identical(attributes(x)$class, "VAR1varyLjBFA")
}