#define ARMA_DONT_PRINT_ERRORS //So the cholesky warning is suppressed
#include <iostream>
#include <fstream>
#include <string>
#include <RcppArmadillo.h>
#include "MCMC_bfaSpatemp.h"

//Function to sample latent polya-gamma process using Gibbs sampling step------------------------------------------------
arma::mat SampleOmega(int f, VAR1datobjFixedL DatObj, VAR1paraFixedL Para) {
  
  //Set data objects
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  int N = DatObj.N;
  int O = DatObj.O;
  arma::mat TrialsMat = DatObj.Trials(arma::span::all, arma::span(f, f), arma::span::all);
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  
  //Moments
  arma::cube MeanOut(N, 1, 1);
  MeanOut(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = Mean;
  MeanOut = arma::reshape(MeanOut, M, O, Nu);
  arma::mat MeanMat = MeanOut(arma::span::all, arma::span(f, f), arma::span::all);

  //Sample latent Variable from full conditional
  arma::mat omega = arma::reshape(pgRcpp(arma::vectorise(TrialsMat), arma::vectorise(MeanMat)), M, Nu);
  return omega;
  
}



//Function to sample latent probit process using Gibbs sampling step------------------------------------------------
VAR1datobjFixedL SampleUpper(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, dataug DatAug) {
  
  //Set data objects
  arma::colvec YStar = DatObj.YStar;
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  arma::cube Cov = Para.Cov;
  
  //Set data augmentation objects
  arma::uvec WhichAbove = DatAug.WhichAbove;
  int NAbove = DatAug.NAbove;
  
  //Moments
  arma::colvec Mu = Mean(WhichAbove);
  arma::colvec SDFull = arma::sqrt(arma::vectorise(Cov));
  arma::colvec SD = SDFull(WhichAbove);
  
  //Sample latent Variable from full conditional
  for (arma::uword i = 0; i < NAbove; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), false, 0);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), 0, arma::datum::inf);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Probit sampling step. Most likely cause for this error is that the data being used is inappropriate (i.e., too far from zero) for a Probit model. Consider scaling towards zero and re-running.");
    YStar(WhichAbove(i)) = Temp;
  }
  DatObj.YStar = YStar;
  return DatObj;
  
}



//Function to sample latent tobit|probit process using Gibbs sampling step------------------------------------------------
VAR1datobjFixedL SampleLower(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, dataug DatAug) {
  
  //Set data objects
  arma::colvec YStar = DatObj.YStar;
  
  //Set parameters
  arma::colvec Mean = Para.Mean;
  arma::cube Cov = Para.Cov;
  
  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::uvec WhichBelow = DatAug.WhichBelow;
  
  //Moments
  arma::colvec Mu = Mean(WhichBelow);
  arma::colvec SDFull = arma::sqrt(arma::vectorise(Cov));
  arma::colvec SD = SDFull(WhichBelow);
  
  //Sample latent Variable from full conditional
  for (arma::uword i = 0; i < NBelow; i++) {
    double Temp = rtnormRcpp(Mu(i), SD(i), true, 0);
    if (!arma::is_finite(Temp)) Temp = rtnormRcppMSM(Mu(i), SD(i), -arma::datum::inf, 0);
    if (!arma::is_finite(Temp)) Rcpp::stop("infinte value sampled in Tobit/Probit sampling step. Most likely cause for this error is that the data being used is inappropriate (i.e., too far from zero) for a Tobit model. Consider scaling towards zero and re-running.");
    YStar(WhichBelow(i)) = Temp;
  }
  DatObj.YStar = YStar;
  return DatObj;
  
}



//Function to sample latent process from its full conditional------------------------------------------------------
std::pair<VAR1datobjFixedL, VAR1paraFixedL> SampleY(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, dataug DatAug) {//called when if (any(FamilyInd != 0))
  
  //Set data objects
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int N = DatObj.N;
  int M = DatObj.M;
  int O = DatObj.O;
  int Nu = DatObj.Nu;
  arma::cube Chi = DatObj.Chi;
  
  //Set parameter objects
  arma::cube Cov = Para.Cov;
  
  //Begin by updating all Tobit/Probit latent variables
  if (DatAug.NBelow > 0) DatObj = SampleLower(DatObj, Para, DatAug); //Tobit|Probit
  if (DatAug.NAbove > 0) DatObj = SampleUpper(DatObj, Para, DatAug); //Probit
  
  if (any(FamilyInd == 3)) {
    arma::colvec YStar = DatObj.YStar;
    //Update Polya-Gamma latent variables
    arma::cube YStarOut(N, 1, 1);
    YStarOut(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = YStar;
    YStarOut = arma::reshape(YStarOut, M, O, Nu);
  
    //Loop over latent dimensions
    for (arma::uword f = 0; f < O; f++) {
  
      //Observation specifics
      int FamInd = FamilyInd(f);
  
      //Polya-Gamma latent process updates
      if (FamInd == 3) {
        
        //Sample omega
        arma::mat chi = Chi(arma::span::all, arma::span(f, f), arma::span::all);
        arma::mat omega = SampleOmega(f, DatObj, Para); //M x Nu
        arma::mat omegaInv = arma::pow(omega, -1);
        
        //Update Polya-Gamma objects
        YStarOut(arma::span::all, arma::span(f, f), arma::span::all) = (chi % omegaInv);
        Cov(arma::span::all, arma::span(f, f), arma::span::all) = omegaInv;
  
      }
      
    //End loop over observation dimension
    }
    YStar = arma::vectorise(YStarOut);
    //Save output
    DatObj.YStar = YStar;
    Para.Cov = Cov;
  }
  
  DatObj.YStarWide = arma::reshape(DatObj.YStar, M * O, Nu);
  return std::pair<VAR1datobjFixedL, VAR1paraFixedL>(DatObj, Para);

}



//Function to sample sigma2(s_i) using a Gibbs sampler step---------------------------------------------------------------
VAR1paraFixedL SampleSigma2(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara) {
  
  //Set data objects
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int Nu = DatObj.Nu;
  int O = DatObj.O;
  int M = DatObj.M;
  int C = DatObj.C;
  arma::mat YStarWide = DatObj.YStarWide;
  
  //Set parameter objects
  arma::colvec Mean = Para.Mean;
  arma::mat MeanMat = arma::reshape(Mean, M * O, Nu);
  arma::mat Sigma2 = Para.Sigma2;
  arma::cube Cov = Para.Cov;

  //Set hyperparameter objects
  double A = HyPara.A;
  double B = HyPara.B;
  
  //Shape is constant over locations
  double Shape = A + 0.5 * Nu;
  
  //Declarations
  arma::uword Index;
  
  //Loop over observations
  arma::uvec NotCount = find(FamilyInd != 3);
  for (arma::uword c = 0; c < (O - C); c++) {
    
    //Set the observation 
    int o = NotCount(c);
    
    //Loop over locations
    for (arma::uword i = 0; i < M; i++) {
    
      //Location specific objects
      Index = i + M * o;
      arma::rowvec Diff = YStarWide.row(Index) - MeanMat.row(Index);
  
      //Calculate rate
      double Resids = arma::as_scalar(Diff * arma::trans(Diff));
    
      //Sample sigma2Io
      double Rate = B + 0.5 * Resids;
      Sigma2(i, c) = rigammaRcpp(Shape, Rate);
      
    //End loop over locations  
    }
    
    //Update covariance
    Cov(arma::span::all, arma::span(o, o), arma::span::all) = arma::repmat(Sigma2.col(c), 1, Nu);

  //End loop over observations 
  }
  
  //Update parameters object
  Para.Sigma2 = Sigma2;
  Para.Cov = Cov;
  return Para;
  
}



//Function to sample beta using a Gibbs sampler step---------------------------------------------------------------
VAR1paraFixedL SampleBeta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara) {
  
  //Set data objects
  arma::mat YStarWide = DatObj.YStarWide;
  int Nu = DatObj.Nu;
  arma::mat EyeNu(Nu, Nu, arma::fill::eye);
  int P = DatObj.P;
  arma::mat X = DatObj.X;
  arma::Col<int> Indeces = DatObj.Indeces;
  
  //Set parameters
  arma::mat Lambda = Para.Lambda;
  arma::cube Cov = Para.Cov;
  arma::colvec Eta = Para.Eta;
  arma::mat BigPhi = Para.BigPhi;
  
  //Set hyperparameters
  arma::mat SigmaBetaInv = HyPara.SigmaBetaInv;
  arma::colvec SigmaBetaInvMuBeta = HyPara.SigmaBetaInvMuBeta;

  //Compute moments
  arma::mat Sum1(P, P, arma::fill::zeros);
  arma::colvec Sum2(P, arma::fill::zeros);
  for (arma::uword t = 0; t < Nu; t++) {
    arma::mat SigmaTInv = arma::diagmat(arma::vectorise(1 / Cov.slice(t)));
    arma::mat XT = X.rows(find(Indeces == t));
    arma::mat tXTSigmaTInv = arma::trans(XT) * SigmaTInv;
    Sum1 += tXTSigmaTInv * XT;
    Sum2 += tXTSigmaTInv * (YStarWide.col(t) - Lambda * BigPhi.col(t));
  }
  
  //Sample Beta
  arma::mat CovBeta = CholInv(Sum1 + SigmaBetaInv);
  arma::colvec MeanBeta = CovBeta * (Sum2 + SigmaBetaInvMuBeta);
  arma::colvec Beta = rmvnormRcpp(1, MeanBeta, CovBeta);
  
  //Update parameters dependent on delta
  arma::colvec XBeta = X * Beta;
  arma::colvec Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;

  //Update parameters object
  Para.Beta = Beta;
  Para.Mean = Mean;
  Para.XBeta = XBeta;
  return Para;
  
}



//Function to sample new value of A using a Metropolis sampler step-----------------------------------------------
VAR1paraFixedL SampleA(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara) {
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    arma::mat eyeK(K, K, arma::fill::eye);
    arma::mat Upsilon = Para.Upsilon;
    arma::mat BigPhi = Para.BigPhi; // k x T
    arma::mat BigPhishedcolT = BigPhi;
    arma::mat V = HyPara.V; //k x k
    arma::colvec Mvec = HyPara.Mvec;//of length k^2
    arma::mat Vinv = CholInv(V);
    BigPhishedcolT.shed_col(Nu - 1);
    arma::mat VA = CholInv(BigPhishedcolT * arma::trans(BigPhishedcolT) + Vinv); // note that \eta_0 is a vector of length k
    arma::colvec meanVec = arma::kron(VA * Vinv, eyeK) * Mvec;
    arma::colvec etaTminus1(K, arma::fill::zeros), etaT(K);
    for (arma::uword t = 0; t < Nu; t++) {
        if (t > 0) { etaTminus1 = etaT; }
        etaT = BigPhi.col(t);
        meanVec += arma::kron(VA * etaTminus1, eyeK) * etaT;
    }
    arma::mat CovMat = arma::kron(VA, Upsilon);
    arma::colvec avec = rmvnormRcpp(1, meanVec, CovMat);
    arma::mat A = arma::reshape(avec, K, K);//k x k
    Para.A = A;
    return Para;
}



//Function to sample Upsilon using a Gibbs sampler step-------------------------------------------------------------------
VAR1paraFixedL SampleUpsilon(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara) {
  
  //Set data objects
  int Nu = DatObj.Nu;
  int K = DatObj.K;

  //Set parameters
  arma::mat BigPhi = Para.BigPhi; //k x T
  arma::mat A = Para.A;

  //Set hyperparameter objects
  double Zeta = HyPara.Zeta;
  arma::mat Omega = HyPara.Omega;
  
  //Compute compSum
  arma::mat compSum(K, K, arma::fill::zeros);
  arma::colvec etaT = BigPhi.col(0); //eta1
  compSum += etaT * arma::trans(etaT);
  arma::colvec etaTminus1(K);
  for (arma::uword t = 1; t < Nu; t++){
      etaTminus1 = etaT;
      etaT = BigPhi.col(t);
      arma::colvec summandVec = etaT - A * etaTminus1;
      compSum += summandVec * arma::trans(summandVec);
  }
  
  //Sample Upsilon
  double n = Zeta + Nu;
  arma::mat V = compSum + Omega;
  arma::mat Upsilon(K, K), UpsilonInv(K, K);
  if (K > 1) {
    UpsilonInv = rwishRcpp(n, CholInv(V));
    Upsilon = CholInv(UpsilonInv);
  } 
  else {
    Upsilon(0, 0) = rigammaRcpp(0.5 * n, 0.5 * arma::as_scalar(V));
    UpsilonInv = 1 / Upsilon;
  }
  
  //Update parameters object
  Para.Upsilon = Upsilon;
  Para.UpsilonInv = UpsilonInv;
  return Para;
}



//Function to sample eta using a Gibbs sampler step---------------------------------------------------------------
VAR1paraFixedL SampleEta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara) {
  //Set data objects
  int K = DatObj.K;
  arma::mat EyeK(K, K, arma::fill::eye);
  arma::mat YStarWide = DatObj.YStarWide; 
  int Nu = DatObj.Nu;
  arma::mat EyeNu(Nu, Nu, arma::fill::eye);
  int M = DatObj.M;
  int O = DatObj.O;
      
  //Set parameters
  arma::mat BigPhi = Para.BigPhi;
  arma::mat Lambda = Para.Lambda;
  arma::cube Cov = Para.Cov;
  arma::mat UpsilonInv = Para.UpsilonInv;
  arma::mat A = Para.A;
  arma::colvec XBeta = Para.XBeta;
  arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);
  
  //Sample etaT
  arma::mat SigmaInv = arma::diagmat(arma::vectorise(1 / Cov.slice(0)));
  arma::mat tLambdaSigmaInv = arma::trans(Lambda) * SigmaInv;
  arma::mat CovEta = CholInv(tLambdaSigmaInv * Lambda + UpsilonInv + arma::trans(A) * UpsilonInv * A);
  arma::mat cholSigma;
  try {
      cholSigma = arma::chol(CovEta);
  }
  catch (...) {
      cholSigma = getCholRobust(CovEta);
  }
  arma::colvec etaTminus1(K, arma::fill::zeros), etaTplus1(K), etaT(K), MeanEtaT(K);
  for (arma::uword t = 0; t < (Nu -  1); t++) { //for time = 1, 2, ..., (T - 1)
      if (t > 0) {etaTminus1 = BigPhi.col(t - 1);}
      etaTplus1 = BigPhi.col(t + 1);
      MeanEtaT = CovEta * (tLambdaSigmaInv * (YStarWide.col(t) - XBetaMat.col(t)) + UpsilonInv * A * etaTminus1 + arma::trans(A) * UpsilonInv * etaTplus1);
      etaT = rmvnormRcppNew(1, MeanEtaT, cholSigma);
      BigPhi.col(t) = etaT;
      //End loop over visits  
  }
  //for time = T
  etaTminus1 = BigPhi.col(Nu - 2);
  CovEta = CholInv(tLambdaSigmaInv * Lambda + UpsilonInv);
  MeanEtaT = CovEta * (tLambdaSigmaInv * (YStarWide.col(Nu - 1) - XBetaMat.col(Nu - 1)) + UpsilonInv * A * etaTminus1);
  etaT = rmvnormRcpp(1, MeanEtaT, CovEta);
  BigPhi.col(Nu - 1) = etaT;
                        
  //Update parameters dependent on eta
  arma::colvec Eta = arma::vectorise(BigPhi);  
  //Update parameters object
  Para.Eta = Eta;
  Para.BigPhi = BigPhi;
  Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
  return Para;
}



//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
paraCLfixedL SampleDelta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara, paraCLfixedL ParaCL) {
    //Set data objects
    int K = DatObj.K;
    int L = DatObj.L;
    int GS = DatObj.GS;
    //Set parameter objects
    arma::mat Theta = ParaCL.Theta;
    arma::colvec Delta = ParaCL.Delta;
    //Set hyperparameter objects
    double A1 = HyPara.A1;
    double A2 = HyPara.A2;
    //Gamma shrinkage prior
    if (GS == 1) {
        double AH;
        for (arma::uword h = 0; h < K; h++) {
            arma::colvec DeltaMinusH = Delta;
            DeltaMinusH(h) = 1;
            //Asign hyperparameter
            if (h == 0) AH = A1;
            else if (h > 0) AH = A2;
            double Resids = 0;
            for (arma::uword j = h; j < K; j++) {
                arma::colvec ThetaJ = Theta.col(j);
                double tThetaTheta = arma::as_scalar(arma::trans(ThetaJ) * ThetaJ);
                Resids += tThetaTheta * arma::prod(DeltaMinusH(arma::span(0, j)));
            }
            double Shape = AH;
            Shape += 0.5 * (K - h) * L;
            double Rate = 1 + 0.5 * Resids;
            Delta(h) = rgammaRcpp(Shape, Rate);
        }
    }
    //NO Gamma shrinkage prior
    else {
        for (arma::uword j = 0; j < K; j++) {
            arma::colvec ThetaJ = Theta.col(j);
            double tThetaTheta = arma::as_scalar(arma::trans(ThetaJ) * ThetaJ);
            double Shape = A1 + 0.5 * L;
            double Rate = A2 + 0.5 * tThetaTheta;
            Delta(j) = rgammaRcpp(Shape, Rate); 
        }
    }
    //Update parameters object
    ParaCL.Delta = Delta;
    if (GS == 1) { ParaCL.Tau = arma::cumprod(Delta); }
    else { ParaCL.Tau = Delta; } //if (GS == 0) 
    return ParaCL;  
}



//Functions to sample new value of rho using a Metropolis sampler step-----------------------------------------------
std::pair<spatpara1, VAR1metrobj> SampleRho(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj) {
    int SpCorInd = DatObj.SpCorInd;
    if (SpCorInd == 0) {//spatial.structure = "continuous"
        //Set data objects
        int M = DatObj.M;
        int O = DatObj.O;
        arma::mat SpDist = DatObj.SpDist;
        int L = DatObj.L;
        int K = DatObj.K;
        arma::colvec ZeroOM(O * M, arma::fill::zeros);
        //Set parameter objects
        double Rho = SpatPara.Rho;
        arma::mat RootiKappa = Para.RootiKappa;
        arma::mat RootiSpCov = SpatPara.RootiSpCov;

        //Set hyperparameter objects
        double ARho = HyPara.ARho;
        double BRho = HyPara.BRho;
        //Set metropolis objects
        double MetropRho = sqrt(MetrObj.MetropRho);
        double AcceptanceRho = MetrObj.AcceptanceRho;
        //Transform current state to real line
        double BigDelta = log((Rho - ARho) / (BRho - Rho));

        //Numerical fix for when the propopsal cholesky doesn't exist
        double RhoProposal, BigDeltaProposal;
        arma::mat CholSpCovProposal(M, M), SpCovProposal(M, M);
        bool Cholesky = false;
        while (!Cholesky) {
            //Sample a new Proposal
            BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropRho));
            //Compute Phi Proposal
            RhoProposal = (BRho * exp(BigDeltaProposal) + ARho) / (1 + exp(BigDeltaProposal));
            //Proposal spatial correlation
            SpCovProposal = SpEXP(RhoProposal, SpDist);
            Cholesky = arma::chol(CholSpCovProposal, SpCovProposal);
        }
        arma::mat RootiSpCovProposal = arma::inv(arma::trimatu(CholSpCovProposal));

        //Alpha structure components
        double Component1A = 0, Component1B = 0;
        arma::mat RootiAlpha, RootiAlphaProposal;
        if (DatObj.CL == 1) {
            RootiAlpha = arma::kron(RootiSpCov, RootiKappa);
            RootiAlphaProposal = arma::kron(RootiSpCovProposal, RootiKappa);
            arma::cube Alpha = ParaCL.Alpha;
            //Loop over columns of j
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                //Loop over mixture components
                for (arma::uword l = 0; l < L - 1; l++) {
                    arma::colvec AlphaJL = AlphaJ.row(l).t();
                    Component1A += lndMvn(AlphaJL, ZeroOM, RootiAlphaProposal);
                    Component1B += lndMvn(AlphaJL, ZeroOM, RootiAlpha);
                }
            }
        }
        else {//if DatObj.CL == 0 
            RootiAlpha = arma::kron(RootiKappa, RootiSpCov);
            RootiAlphaProposal = arma::kron(RootiKappa, RootiSpCovProposal);
            arma::mat Lambda = Para.Lambda;
            //Loop over columns of j
            for (arma::uword j = 0; j < K; j++) {
                arma::colvec LambdaJ = Lambda.col(j);
                Component1A += lndMvn(LambdaJ, ZeroOM, RootiAlphaProposal);
                Component1B += lndMvn(LambdaJ, ZeroOM, RootiAlpha);                
            }
        }

        double Component1 = Component1A - Component1B;
        double Component2A = BigDeltaProposal;
        double Component2B = BigDelta;
        double Component2 = Component2A - Component2B;
        double Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
        //Compute log acceptance ratio
        double LogR = Component1 + Component2 + Component3;

        //Metropolis update
        double RandU = randuRcpp();
        if (log(RandU) < LogR) {
            //Keep Count of Acceptances
            AcceptanceRho++;
            MetrObj.AcceptanceRho = AcceptanceRho;
            //Update parameters object
            SpatPara.Rho = RhoProposal;
            SpatPara.RootiSpCov = RootiSpCovProposal;
            SpatPara.SpCovInv = RootiSpCovProposal * arma::trans(RootiSpCovProposal);
        }
    }
    else { //when SpCorInd == 1, i.e., spatial.structure = "discrete"
        //Set data objects
        int M = DatObj.M;
        int O = DatObj.O;
        arma::mat SpDist = DatObj.SpDist;
        int L = DatObj.L;
        int K = DatObj.K;
        arma::colvec ZeroOM(O * M, arma::fill::zeros); 
        //Set parameter objects
        double Rho = SpatPara.Rho;
        arma::mat KappaInv = Para.KappaInv;
        arma::mat RootiKappa = arma::trans(getCholRobust(KappaInv)); //transpose of the cholesky factor of KappaInv instead; not the real rooti but corresponds to RootiSpCov here
        arma::mat RootiSpCov = SpatPara.RootiSpCov; //transpose of the cholesky factor of SpCovInv instead; not the real rooti but works here too
        arma::colvec DwVec = SpatPara.DwVec;
        arma::mat Dw = arma::diagmat(DwVec);

        //Set hyperparameter objects
        double ARho = HyPara.ARho;
        double BRho = HyPara.BRho;
        //Set metropolis objects
        double MetropRho = sqrt(MetrObj.MetropRho);
        double AcceptanceRho = MetrObj.AcceptanceRho;
        //Transform current state to real line
        double BigDelta = log((Rho - ARho) / (BRho - Rho));

        //Numerical fix for when the propopsal cholesky doesn't exist
        double RhoProposal, BigDeltaProposal;
        arma::mat RootiSpCovProposal(M, M), SpCovInvProposal(M, M);
        bool Cholesky = false;
        while (!Cholesky) {
            //Sample a new Proposal
            BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropRho));
            //Compute Phi Proposal
            RhoProposal = (BRho * exp(BigDeltaProposal) + ARho) / (1 + exp(BigDeltaProposal));
            //Proposal spatial correlation
            SpCovInvProposal = Dw - RhoProposal * SpDist;
            Cholesky = arma::chol(RootiSpCovProposal, SpCovInvProposal);
        }
        RootiSpCovProposal = arma::trans(RootiSpCovProposal);

        //Alpha structure components
        double Component1A = 0, Component1B = 0;
        arma::mat RootiAlpha, RootiAlphaProposal;
        if (DatObj.CL == 1) {
            RootiAlpha = arma::kron(RootiSpCov, RootiKappa);
            RootiAlphaProposal = arma::kron(RootiSpCovProposal, RootiKappa);
            arma::cube Alpha = ParaCL.Alpha;
            //Loop over columns of j
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                //Loop over mixture components
                for (arma::uword l = 0; l < L - 1; l++) {
                    arma::colvec AlphaJL = AlphaJ.row(l).t();
                    Component1A += lndMvn(AlphaJL, ZeroOM, RootiAlphaProposal);
                    Component1B += lndMvn(AlphaJL, ZeroOM, RootiAlpha);
                }
            }
        }
        else {//if DatObj.CL == 0 
            RootiAlpha = arma::kron(RootiKappa, RootiSpCov);
            RootiAlphaProposal = arma::kron(RootiKappa, RootiSpCovProposal);
            arma::mat Lambda = Para.Lambda;
            //Loop over columns of j
            for (arma::uword j = 0; j < K; j++) {
                arma::colvec LambdaJ = Lambda.col(j);
                Component1A += lndMvn(LambdaJ, ZeroOM, RootiAlphaProposal);
                Component1B += lndMvn(LambdaJ, ZeroOM, RootiAlpha);
            }
        }

        double Component1 = Component1A - Component1B;
        double Component2A = BigDeltaProposal;
        double Component2B = BigDelta;
        double Component2 = Component2A - Component2B;
        double Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
        //Compute log acceptance ratio
        double LogR = Component1 + Component2 + Component3;

        //Metropolis update
        double RandU = randuRcpp();
        if (log(RandU) < LogR) {
            //Keep Count of Acceptances
            AcceptanceRho++;
            MetrObj.AcceptanceRho = AcceptanceRho;
            //Update parameters object
            SpatPara.Rho = RhoProposal;
            SpatPara.RootiSpCov = RootiSpCovProposal;
            SpatPara.SpCovInv = SpCovInvProposal;
        }
    }
    //Return output object
    return std::pair<spatpara1, VAR1metrobj>(SpatPara, MetrObj);
}

std::pair<spatpara2, VAR1metrobj> SampleRho(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj) {//when IS == 1 && SpCorInd == 0 && DatObj.SA == 1
    //Set data objects
    int M = DatObj.M;
    int O = DatObj.O;
    arma::mat SpDist = DatObj.SpDist;
    int L = DatObj.L;
    int K = DatObj.K;
    int h = DatObj.h;
    //Set parameter objects
    double Rho = SpatPara.Rho;
    arma::mat KappaInv = Para.KappaInv;
    arma::mat SpCovInvApprox = SpatPara.SpCovInv;
    double detSpCovApprox = SpatPara.detSpCov;
    arma::umat nnInd = SpatPara.nnInd;

    //Set hyperparameter objects
    double ARho = HyPara.ARho;
    double BRho = HyPara.BRho;
    //Set metropolis objects
    double MetropRho = sqrt(MetrObj.MetropRho);
    double AcceptanceRho = MetrObj.AcceptanceRho;
    //Transform current state to real line
    double BigDelta = log((Rho - ARho) / (BRho - Rho));

    //Numerical fix for when the propopsal cholesky doesn't exist
    double RhoProposal, BigDeltaProposal;
    double detSpCovApproxProposal = 1;
    arma::mat SpCovInvApproxProposal(M, M, arma::fill::zeros);
    arma::mat SpCovProposal(M, M);
    //Sample a new Proposal
    BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropRho));
    //Compute Phi Proposal
    RhoProposal = (BRho * exp(BigDeltaProposal) + ARho) / (1 + exp(BigDeltaProposal));
    //Proposal spatial correlation
    SpCovProposal = SpEXP(RhoProposal, SpDist);

    //calculate SpCovInvApproxProposal and detSpCovApproxProposal
    int nnMaxNum;
    SpCovInvApproxProposal(0, 0) = 1; // 1 / Fs(0) as BsStarRow0's first entry is 1 and all other entries are 0
    for (int i = 1; i < M; i++) {
        if (i < h) nnMaxNum = i - 1;
        else nnMaxNum = h - 1;
        //arma::Row<arma::uword> nni = nnInd(arma::span(i, i), arma::span(0, nnMaxNum));
        arma::uvec nni = nnInd(arma::span(i, i), arma::span(0, nnMaxNum)).t();
        arma::rowvec SpCovProposalRowi = SpCovProposal.row(i);
        arma::mat SpCovProposalinni = SpCovProposalRowi(nni);
        arma::mat Bsi = SpCovProposalinni.t() * CholInv(SpCovProposal(nni, nni));
        double Fsi = arma::as_scalar(1 - Bsi * SpCovProposalinni);
        detSpCovApproxProposal *= Fsi;
        //nni.insert_cols(nnMaxNum + 1, arma::uvec({ arma::uword(i) }));
        //Bsi.insert_cols(nnMaxNum + 1, arma::vec({ -1 }));
        //SpCovInvApproxProposal(nni, nni) += arma::trans(Bsi) * (1 / Fsi) * (Bsi);
        SpCovInvApproxProposal(nni, nni) += arma::trans(-Bsi) * (1 / Fsi) * (-Bsi);
        arma::rowvec SpCovInvApproxProposalRowi(M, arma::fill::zeros);
        SpCovInvApproxProposalRowi(nni) = -Bsi / Fsi;
        SpCovInvApproxProposalRowi(i) = 1 / Fsi;
        SpCovInvApproxProposal.row(i) = SpCovInvApproxProposalRowi;
        SpCovInvApproxProposal.col(i) = SpCovInvApproxProposalRowi.t();
    }

    //Alpha structure components
    double component1 = 0;
    double logDetDiff = O * log(detSpCovApproxProposal / detSpCovApprox);
    arma::mat SpCovInvDiff = SpCovInvApproxProposal - SpCovInvApprox;
    //arma::mat kronInvDiff = arma::kron(SpCovInvDiff, KappaInv);
    //arma::mat kronInvDiff2 = arma::kron(KappaInv, SpCovInvDiff);
    arma::sp_mat kronInvDiff(arma::kron(SpCovInvDiff, KappaInv));
    arma::sp_mat kronInvDiff2(arma::kron(KappaInv, SpCovInvDiff));

    if (DatObj.CL == 1) {//alphaMethod must = "block" here if DatObj.CL == 1
        component1 -= logDetDiff * (L - 1) * K;
        arma::cube Alpha = ParaCL.Alpha;
        //Loop over columns of j
        for (arma::uword j = 0; j < K; j++) {           
            arma::mat AlphaJ = Alpha.slice(j);
            //Loop over mixture components
            for (arma::uword l = 0; l < L - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                component1 -= arma::as_scalar(AlphaJL * kronInvDiff * arma::trans(AlphaJL));
            }
        }
        component1 /= 2;
    }
    else {//if DatObj.CL == 0 so that L = 1
        component1 -= logDetDiff * K;
        arma::mat Lambda = Para.Lambda;
        //Loop over columns of j
        for (arma::uword j = 0; j < K; j++) {
            //Objects that depend on j           
            arma::colvec LambdaJ = Lambda.col(j);
            component1 -= arma::as_scalar(arma::trans(LambdaJ) * kronInvDiff2 * LambdaJ);
        }
        component1 /= 2;
    }

    double Component2A = BigDeltaProposal;
    double Component2B = BigDelta;
    double Component2 = Component2A - Component2B;
    double Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
    //Compute log acceptance ratio
    double LogR = component1 + Component2 + Component3;

    //Metropolis update
    double RandU = randuRcpp();
    if (log(RandU) < LogR) {
        //Keep Count of Acceptances
        AcceptanceRho++;
        MetrObj.AcceptanceRho = AcceptanceRho;
        //Update parameters object
        SpatPara.Rho = RhoProposal;
        SpatPara.detSpCov = detSpCovApproxProposal;
        SpatPara.SpCovInv = SpCovInvApproxProposal;
    }

    //Return output object
    return std::pair<spatpara2, VAR1metrobj>(SpatPara, MetrObj);
}

std::pair<spatpara3, VAR1metrobj> SampleRho(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj) {//when IS == 1 && SpCorInd == 0 && DatObj.SA == 1, DatObj.CL == 1, alphaMethod == "sequential"
        //Set data objects
        int M = DatObj.M;
        int O = DatObj.O;
        arma::mat SpDist = DatObj.SpDist;
        int L = DatObj.L;
        int K = DatObj.K;
        int h = DatObj.h;
        //Set parameter objects
        double Rho = SpatPara.Rho;
        arma::mat KappaInv = Para.KappaInv;
        arma::mat SpCovInvApprox = SpatPara.SpCovInv;
        double detSpCovApprox = SpatPara.detSpCov;
        arma::cube Alpha = ParaCL.Alpha;
        arma::umat nnInd = SpatPara.nnInd;

        //Set hyperparameter objects
        double ARho = HyPara.ARho;
        double BRho = HyPara.BRho;
        //Set metropolis objects
        double MetropRho = sqrt(MetrObj.MetropRho);
        double AcceptanceRho = MetrObj.AcceptanceRho;
        //Transform current state to real line
        double BigDelta = log((Rho - ARho) / (BRho - Rho));

        //Numerical fix for when the propopsal cholesky doesn't exist
        double RhoProposal, BigDeltaProposal;
        double detSpCovApproxProposal = 1;
        arma::mat SpCovInvApproxProposal(M, M, arma::fill::zeros);
        arma::mat SpCovProposal(M, M);
        //Sample a new Proposal
        BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropRho));
        //Compute Phi Proposal
        RhoProposal = (BRho * exp(BigDeltaProposal) + ARho) / (1 + exp(BigDeltaProposal));
        //Proposal spatial correlation
        SpCovProposal = SpEXP(RhoProposal, SpDist);

        //calculate SpCovInvApproxProposal, detSpCovApproxProposal, and Fs, Bs
        arma::colvec Fs(M, arma::fill::ones);
        arma::mat Bs(M, h, arma::fill::zeros);
        SpCovInvApproxProposal(0, 0) = 1; // 1 / Fs(0) as BsStarRow0's first entry is 1 and all other entries are 0
        int nnMaxNum;
        for (int i = 1; i < M; i++) {
            if (i < h) nnMaxNum = i - 1;
            else nnMaxNum = h - 1;
            arma::uvec nni = nnInd(arma::span(i, i), arma::span(0, nnMaxNum)).t();
            arma::rowvec SpCovProposalRowi = SpCovProposal.row(i);
            arma::mat SpCovProposalinni = SpCovProposalRowi(nni);
            arma::mat Bsi = SpCovProposalinni.t() * CholInv(SpCovProposal(nni, nni));
            Bs(arma::span(i, i), arma::span(0, nnMaxNum)) = Bsi;
            double Fsi = arma::as_scalar(1 - Bsi * SpCovProposalinni);
            Fs(i) = Fsi;
            detSpCovApproxProposal *= Fsi;
            //nni.insert_cols(nnMaxNum + 1, arma::uvec({ arma::uword(i) }));
            //Bsi.insert_cols(nnMaxNum + 1, arma::vec({ -1 }));
            //SpCovInvApproxProposal(nni, nni) += arma::trans(Bsi) * (1 / Fsi) * (Bsi);
            SpCovInvApproxProposal(nni, nni) += arma::trans(-Bsi) * (1 / Fsi) * (-Bsi);
            arma::rowvec SpCovInvApproxProposalRowi(M, arma::fill::zeros);
            SpCovInvApproxProposalRowi(nni) = -Bsi / Fsi;
            SpCovInvApproxProposalRowi(i) = 1 / Fsi;
            SpCovInvApproxProposal.row(i) = SpCovInvApproxProposalRowi;
            SpCovInvApproxProposal.col(i) = SpCovInvApproxProposalRowi.t();
        }

        //Alpha structure components
        double component1 = 0;
        double logDetDiff = O * log(detSpCovApproxProposal / detSpCovApprox);
        arma::sp_mat kronInvDiff(arma::kron((SpCovInvApproxProposal - SpCovInvApprox), KappaInv));
        //arma::mat kronInvDiff = arma::kron((SpCovInvApproxProposal - SpCovInvApprox), KappaInv);
        
        //Loop over columns of j (DatObj.CL must = 1 here)
        for (arma::uword j = 0; j < K; j++) {
            component1 -= logDetDiff * (L - 1);
            arma::mat AlphaJ = Alpha.slice(j);
            //Loop over mixture components
            for (arma::uword l = 0; l < L - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                component1 -= arma::as_scalar(AlphaJL * kronInvDiff * arma::trans(AlphaJL));
            }
        }
        component1 /= 2;

        double Component2A = BigDeltaProposal;
        double Component2B = BigDelta;
        double Component2 = Component2A - Component2B;
        double Component3 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
        //Compute log acceptance ratio
        double LogR = component1 + Component2 + Component3;

        //Metropolis update
        double RandU = randuRcpp();
        if (log(RandU) < LogR) {
            //Keep Count of Acceptances
            AcceptanceRho++;
            MetrObj.AcceptanceRho = AcceptanceRho;
            //Update parameters object
            SpatPara.Rho = RhoProposal;
            SpatPara.detSpCov = detSpCovApproxProposal;
            SpatPara.SpCovInv = SpCovInvApproxProposal;
            SpatPara.Bs = Bs;
            SpatPara.Fs = Fs;
        }
        //Return output object
        return std::pair<spatpara3, VAR1metrobj>(SpatPara, MetrObj);
}



//Functions to sample kappa using a Gibbs sampler step---------------------------------------------------------------
VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, VAR1hypara HyPara) {//CL = 1
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    arma::cube Alpha = ParaCL.Alpha;
    arma::mat SpCovInv = SpatPara.SpCovInv;
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        for (arma::uword l = 0; l < L - 1; l++) {
            arma::rowvec AlphaJL = Alpha.slice(j).row(l);
            arma::mat AlphaJLMat = arma::reshape(AlphaJL, O, M);
            Resids += AlphaJLMat * SpCovInv * arma::trans(AlphaJLMat);
        }
    }
    double Shape = SmallUpsilon + M * (L - 1) * K;
    arma::mat Rate = Resids + BigTheta;
    arma::mat Kappa(O, O), KappaInv(O, O), RootiKappa(O, O);
    if (O == 1) {
        Kappa(0, 0) = rigammaRcpp(0.5 * arma::as_scalar(Shape), 0.5 * arma::as_scalar(Rate));
        KappaInv = 1 / Kappa;
        RootiKappa = arma::sqrt(KappaInv);
    }
    else if (O > 1) {
        KappaInv = rwishRcpp(Shape, CholInv(Rate));
        Kappa = CholInv(KappaInv);
        RootiKappa = GetRooti(Kappa);
    }

    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;
}

VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara1 SpatPara, VAR1hypara HyPara) {//CL = 0
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    arma::mat Lambda = Para.Lambda;
    arma::mat SpCovInv = SpatPara.SpCovInv;
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        arma::colvec LambdaJ = Lambda.col(j);
        arma::mat LambdaJMat = arma::reshape(LambdaJ, M, O);
        Resids += arma::trans(LambdaJMat) * SpCovInv * LambdaJMat;
    }
    double Shape = SmallUpsilon + M * K;
    arma::mat Rate = Resids + BigTheta;
    arma::mat Kappa(O, O), KappaInv(O, O), RootiKappa(O, O);
    if (O == 1) {
        Kappa(0, 0) = rigammaRcpp(0.5 * arma::as_scalar(Shape), 0.5 * arma::as_scalar(Rate));
        KappaInv = 1 / Kappa;
        RootiKappa = arma::sqrt(KappaInv);
    }
    else if (O > 1) {
        KappaInv = rwishRcpp(Shape, CholInv(Rate));
        Kappa = CholInv(KappaInv);
        RootiKappa = GetRooti(Kappa);
    }

    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;
}

VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, VAR1hypara HyPara) {//CL = 1
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    arma::cube Alpha = ParaCL.Alpha;
    arma::sp_mat SpCovInv(SpatPara.SpCovInv);
    //arma::mat SpCovInv = SpatPara.SpCovInv;
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        for (arma::uword l = 0; l < L - 1; l++) {
            arma::rowvec AlphaJL = Alpha.slice(j).row(l);
            arma::mat AlphaJLMat = arma::reshape(AlphaJL, O, M);
            Resids += AlphaJLMat * SpCovInv * arma::trans(AlphaJLMat);        
        }
    }
    double Shape = SmallUpsilon + M * (L - 1) * K;
    arma::mat Rate = Resids + BigTheta;
    arma::mat Kappa(O, O), KappaInv(O, O), RootiKappa(O, O);
    if (O == 1) {
        Kappa(0, 0) = rigammaRcpp(0.5 * arma::as_scalar(Shape), 0.5 * arma::as_scalar(Rate));
        KappaInv = 1 / Kappa;
        RootiKappa = arma::sqrt(KappaInv);
    }
    else if (O > 1) {
        KappaInv = rwishRcpp(Shape, CholInv(Rate));
        Kappa = CholInv(KappaInv);
        RootiKappa = GetRooti(Kappa);
    }

    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;
}

VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara2 SpatPara, VAR1hypara HyPara) {//CL = 0
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    arma::mat Lambda = Para.Lambda;
    arma::sp_mat SpCovInv(SpatPara.SpCovInv);
    //arma::mat SpCovInv = SpatPara.SpCovInv;
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        arma::colvec LambdaJ = Lambda.col(j);
        arma::mat LambdaJMat = arma::reshape(LambdaJ, M, O);
        Resids += arma::trans(LambdaJMat) * SpCovInv * LambdaJMat;    
    }
    double Shape = SmallUpsilon + M * K;
    arma::mat Rate = Resids + BigTheta;
    arma::mat Kappa(O, O), KappaInv(O, O), RootiKappa(O, O);
    if (O == 1) {
        Kappa(0, 0) = rigammaRcpp(0.5 * arma::as_scalar(Shape), 0.5 * arma::as_scalar(Rate));
        KappaInv = 1 / Kappa;
        RootiKappa = arma::sqrt(KappaInv);
    }
    else if (O > 1) {
        KappaInv = rwishRcpp(Shape, CholInv(Rate));
        Kappa = CholInv(KappaInv);
        RootiKappa = GetRooti(Kappa);
    }

    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;
}

VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, VAR1hypara HyPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    arma::cube Alpha = ParaCL.Alpha;
    arma::sp_mat SpCovInv(SpatPara.SpCovInv);
    //arma::mat SpCovInv = SpatPara.SpCovInv;
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        for (arma::uword l = 0; l < L - 1; l++) {
            arma::rowvec AlphaJL = Alpha.slice(j).row(l);
            arma::mat AlphaJLMat = arma::reshape(AlphaJL, O, M);
            Resids += AlphaJLMat * SpCovInv * arma::trans(AlphaJLMat);
        }
    }
    double Shape = SmallUpsilon + M * (L - 1) * K;
    arma::mat Rate = Resids + BigTheta;
    arma::mat Kappa(O, O), KappaInv(O, O), RootiKappa(O, O);
    if (O == 1) {
        Kappa(0, 0) = rigammaRcpp(0.5 * arma::as_scalar(Shape), 0.5 * arma::as_scalar(Rate));
        KappaInv = 1 / Kappa;
        RootiKappa = arma::sqrt(KappaInv);
    }
    else if (O > 1) {
        KappaInv = rwishRcpp(Shape, CholInv(Rate));
        Kappa = CholInv(KappaInv);
        RootiKappa = GetRooti(Kappa);
    }

    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;
}



//Functions to sample alphajl(s_i)'s (when CL = 1) / lambdajio's (when CL = 0) using a Gibbs sampler step---------------------------------------------------------------
paraCLfixedL SampleAlpha(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int keptIter, bool storePostPara) {//CL = 1
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    int alphaWsFiles = DatObj.alphaWsFiles;
    arma::mat EyeOM(O * M, O * M, arma::fill::eye);

    arma::mat KappaInv = Para.KappaInv;
    arma::cube Z = ParaCL.Z;
    arma::cube Alpha = ParaCL.Alpha;
    arma::mat SpCovInv = SpatPara.SpCovInv;

    arma::mat CovAlpha = CholInv(EyeOM + arma::kron(SpCovInv, KappaInv));
    arma::mat cholSigma;
    try {
        cholSigma = arma::chol(CovAlpha);
    }
    catch (...) {
        cholSigma = getCholRobust(CovAlpha);
    }
    for (arma::uword j = 0; j < K; j++) {
        arma::mat AlphaJ(L - 1, M * O);
        arma::mat ZJ = Z.slice(j);
        for (arma::uword l = 0; l < L - 1; l++) {
            arma::colvec MeanAlpha = CovAlpha * arma::trans(ZJ.row(l));
            arma::colvec AlphaJL = rmvnormRcppNew(1, MeanAlpha, cholSigma);
            AlphaJ.row(l) = arma::trans(AlphaJL);
        }
        Alpha.slice(j) = AlphaJ;
    }

    arma::cube logWeights = GetLogWeightsFixedL(Alpha, K, M, L, O);
    ParaCL.Alpha = Alpha;
    ParaCL.logWeights = logWeights;
    if (alphaWsFiles == 1 && storePostPara == true) {
        if (DatObj.spatPred == 1) {
            //std::string filename = "fixedLalphaIter" + to_string(keptIter) + ".txt";
            Alpha.save("fixedLalphaKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
        if (DatObj.storeW == 1) {
            logWeights.save("fixedLlogweightsKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
    }
    if (DatObj.storeW == 1 && storePostPara == true) {
        arma::cube Weights = GetWeightsFixedL(Alpha, K, M, L, O);
        if (alphaWsFiles == 1) {
            Weights.save("fixedLweightsKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }       
        else { // if alphaWsFiles == 0
            ParaCL.Weights = Weights; 
        }
    }
    return ParaCL;        
}

VAR1paraFixedL SampleLambda(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara1 SpatPara) {//CL = 0
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat EyeNu(Nu, Nu, arma::fill::eye);
    arma::mat X = DatObj.X;
    arma::Col<int> Indeces = DatObj.Indeces;
    arma::mat YStarWide = DatObj.YStarWide;

    arma::mat KappaInv = Para.KappaInv;
    arma::mat SpCovInv = SpatPara.SpCovInv;
    arma::mat BigPhi = Para.BigPhi;
    arma::cube Cov = Para.Cov;
    arma::colvec Beta = Para.Beta;
    arma::mat Lambda = Para.Lambda;
    arma::colvec Eta = Para.Eta;
    arma::colvec XBeta = Para.XBeta;

    arma::mat CovAlphaFixed = arma::kron(KappaInv, SpCovInv);//since alpha's here is the columns of Lambda (ordered spatially first and then by observation type)
    arma::mat SigmaTInv = arma::diagmat(arma::vectorise(1 / Cov.slice(0)));
    for (arma::uword j = 0; j < K; j++) {
        arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
        LambdaMinusJ.shed_col(j);
        BigPhiMinusJ.shed_row(j);
        double sum1 = 0;
        arma::colvec sum2(M * O, arma::fill::zeros);
        for (arma::uword t = 0; t < Nu; t++) {            
            arma::mat XT = X.rows(find(Indeces == t));
            double EtaTJ = BigPhi(j, t);
            arma::colvec BigPhiMinusJt = BigPhiMinusJ.col(t);
            arma::colvec LambdaSum = LambdaMinusJ * BigPhiMinusJt;
            arma::colvec MuTJ = (YStarWide.col(t) - XT * Beta - LambdaSum);
            sum1 += std::pow(EtaTJ, 2);
            sum2 += EtaTJ * MuTJ;
        }
        arma::mat CovAlpha = CholInv(sum1 * SigmaTInv + CovAlphaFixed);
        arma::colvec MeanAlpha = CovAlpha * SigmaTInv * sum2;
        arma::colvec AlphaJ = rmvnormRcpp(1, MeanAlpha, CovAlpha);
        Lambda.col(j) = AlphaJ;
    }

    Para.Lambda = Lambda;
    Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
    return Para;
}

paraCLfixedL SampleAlpha(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int keptIter, bool storePostPara) {//CL = 1
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    int alphaWsFiles = DatObj.alphaWsFiles;
    arma::mat EyeOM(O * M, O * M, arma::fill::eye);

    arma::mat KappaInv = Para.KappaInv;
    arma::cube Z = ParaCL.Z;
    arma::cube Alpha = ParaCL.Alpha;
    arma::mat SpCovInv = SpatPara.SpCovInv;

    arma::mat CovAlpha = CholInv(EyeOM + arma::kron(SpCovInv, KappaInv));
    arma::mat cholSigma;
    try {
        cholSigma = arma::chol(CovAlpha);
    }
    catch (...) {
        cholSigma = getCholRobust(CovAlpha);
    }
    for (arma::uword j = 0; j < K; j++) {
        arma::mat AlphaJ(L - 1, M * O);
        arma::mat ZJ = Z.slice(j);
        for (arma::uword l = 0; l < L - 1; l++) {
            arma::colvec MeanAlpha = CovAlpha * arma::trans(ZJ.row(l));
            arma::colvec AlphaJL = rmvnormRcppNew(1, MeanAlpha, cholSigma);
            AlphaJ.row(l) = arma::trans(AlphaJL);
        }
        Alpha.slice(j) = AlphaJ;
    }

    arma::cube logWeights = GetLogWeightsFixedL(Alpha, K, M, L, O);
    ParaCL.Alpha = Alpha;
    ParaCL.logWeights = logWeights;
    if (alphaWsFiles == 1 && storePostPara == true) {
        if (DatObj.spatPred == 1) {
            //std::string filename = "fixedLalphaIter" + to_string(keptIter) + ".txt";
            Alpha.save("fixedLalphaKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
        if (DatObj.storeW == 1) {
            logWeights.save("fixedLlogweightsKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
    }
    if (DatObj.storeW == 1 && storePostPara == true) {
        arma::cube Weights = GetWeightsFixedL(Alpha, K, M, L, O);
        if (alphaWsFiles == 1) {
            Weights.save("fixedLweightsKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
        else { // if alphaWsFiles == 0
            ParaCL.Weights = Weights;
        }
    }
    return ParaCL;
}

VAR1paraFixedL SampleLambda(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara2 SpatPara) {//CL = 0
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat EyeNu(Nu, Nu, arma::fill::eye);
    arma::mat X = DatObj.X;
    arma::Col<int> Indeces = DatObj.Indeces;
    arma::mat YStarWide = DatObj.YStarWide;

    arma::mat KappaInv = Para.KappaInv;
    arma::mat SpCovInv = SpatPara.SpCovInv;
    arma::mat BigPhi = Para.BigPhi;
    arma::cube Cov = Para.Cov;
    arma::colvec Beta = Para.Beta;
    arma::mat Lambda = Para.Lambda;
    arma::colvec Eta = Para.Eta;
    arma::colvec XBeta = Para.XBeta;

    arma::mat CovAlphaFixed = arma::kron(KappaInv, SpCovInv);//since alpha's here is the columns of Lambda (ordered spatially first and then by observation type)
    arma::mat SigmaTInv = arma::diagmat(arma::vectorise(1 / Cov.slice(0)));
    for (arma::uword j = 0; j < K; j++) {
        arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
        LambdaMinusJ.shed_col(j);
        BigPhiMinusJ.shed_row(j);
        double sum1 = 0;
        arma::colvec sum2(M * O, arma::fill::zeros);
        for (arma::uword t = 0; t < Nu; t++) {
            arma::mat XT = X.rows(find(Indeces == t));
            double EtaTJ = BigPhi(j, t);
            arma::colvec BigPhiMinusJt = BigPhiMinusJ.col(t);
            arma::colvec LambdaSum = LambdaMinusJ * BigPhiMinusJt;
            arma::colvec MuTJ = (YStarWide.col(t) - XT * Beta - LambdaSum);
            sum1 += std::pow(EtaTJ, 2);
            sum2 += EtaTJ * MuTJ;
        }
        arma::mat CovAlpha = CholInv(sum1 * SigmaTInv + CovAlphaFixed);
        arma::colvec MeanAlpha = CovAlpha * SigmaTInv * sum2;
        arma::colvec AlphaJ = rmvnormRcpp(1, MeanAlpha, CovAlpha);
        Lambda.col(j) = AlphaJ;
    }

    Para.Lambda = Lambda;
    Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
    return Para;
}

paraCLfixedL SampleAlpha(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int keptIter, bool storePostPara) {
    //DatObj.CL must be 1 in this case (and alphaMethod = "sequential", IS = 1, SpCorInd = 0, spatApprox = TRUE)
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    int h = DatObj.h;
    int alphaWsFiles = DatObj.alphaWsFiles;
    arma::mat EyeO(O, O, arma::fill::eye);
   
    arma::cube Z = ParaCL.Z;
    arma::cube Alpha = ParaCL.Alpha;

    arma::umat nnInd = SpatPara.nnInd;
    arma::mat Bs = SpatPara.Bs;
    arma::colvec Fs = SpatPara.Fs;
    arma::field<arma::mat> whichJs = SpatPara.whichJs;

    if (O == 1) {
        arma::mat KappaInvMat = Para.KappaInv;
        double KappaInv = KappaInvMat(0, 0);
        arma::vec Vs(M);
        for (int i = 0; i < M; i++) {
            double Vsi = 1 + KappaInv / Fs(i);
            if (i < M - 1) {
                arma::mat iAsNN = whichJs(i, 0);
                int iNumJ = size(iAsNN, 1);
                if (!(iNumJ == 1 && iAsNN(0, 0) == 0)) {
                    for (arma::uword nj = 0; nj < iNumJ; nj++) {
                        arma::colvec ijInfo = iAsNN.col(nj);
                        int j = ijInfo(0);
                        int ksi = ijInfo(1);
                        double Bsjksi = Bs(j, ksi);
                        Vsi += std::pow(Bsjksi, 2) / Fs(j) * KappaInv;
                    }
                }
            }
            Vs(i) = 1 / Vsi;
        }
        int nnMaxNum;
        int nnMaxNumJ;
        for (arma::uword jj = 0; jj < K; jj++) {
            arma::mat AlphaJ = Alpha.slice(jj);
            arma::mat ZJ = Z.slice(jj);
            for (arma::uword lj = 0; lj < L - 1; lj++) {
                arma::colvec AlphaJlj = arma::trans(AlphaJ.row(lj));// of length M
                arma::colvec ZJlj = arma::trans(ZJ.row(lj));// of length M
                for (int i = 0; i < M; i++) {
                    double muJljsi = ZJlj(i);
                    if (i > 0) {
                        if (i < h) nnMaxNum = i - 1;
                        else nnMaxNum = h - 1;
                        arma::mat Bsi = Bs(arma::span(i, i), arma::span(0, nnMaxNum));
                        arma::uvec nni = nnInd(arma::span(i, i), arma::span(0, nnMaxNum)).t();
                        arma::colvec AlphaJljNsi = AlphaJlj(nni);
                        muJljsi += KappaInv / Fs(i) * arma::as_scalar(Bsi * AlphaJljNsi);
                    }
                    if (i < M - 1) {
                        arma::mat iAsNN = whichJs(i, 0);
                        int iNumJ = size(iAsNN, 1);
                        if (!(iNumJ == 1 && iAsNN(0, 0) == 0)) {
                            for (arma::uword nj = 0; nj < iNumJ; nj++) {
                                arma::colvec ijInfo = iAsNN.col(nj);
                                int j = ijInfo(0);
                                if (j < h) nnMaxNumJ = j - 1;
                                else nnMaxNumJ = h - 1;
                                int ksi = ijInfo(1);
                                double AlphaJljsj = AlphaJlj(j);
                                arma::rowvec Bsj = Bs(arma::span(j, j), arma::span(0, nnMaxNumJ));
                                double Bsjksi = Bsj(ksi);
                                arma::uvec nnj = nnInd(arma::span(j, j), arma::span(0, nnMaxNumJ)).t();
                                arma::colvec AlphaJljNsj = AlphaJlj(nnj);
                                AlphaJljNsj(ksi) = 0;
                                Bsj(ksi) = 0;
                                muJljsi += (Bsjksi / Fs(j) * KappaInv) * (AlphaJljsj - arma::as_scalar(Bsj * AlphaJljNsj));
                            }
                        }
                    }
                    double Vsi = Vs(i);
                    double AlphaJljsiNew = arma::as_scalar(rnormRcpp(1, Vsi * muJljsi, sqrt(Vsi)));
                    AlphaJlj(i) = AlphaJljsiNew;
                }
                AlphaJ.row(lj) = arma::trans(AlphaJlj); 
            }
            Alpha.slice(jj) = AlphaJ;
        }
    }
    else {
        arma::mat KappaInv = Para.KappaInv;
        arma::cube Vs(O, O, M);
        arma::cube cholVs(O, O, M);
        for (int i = 0; i < M; i++) {
            arma::mat Vsi(O, O, arma::fill::eye);
            Vsi += (1 / Fs(i)) * KappaInv;
            if (i < M - 1) {
                arma::mat iAsNN = whichJs(i, 0);
                int iNumJ = size(iAsNN, 1);
                if (!(iNumJ == 1 && iAsNN(0, 0) == 0)) {
                    for (arma::uword nj = 0; nj < iNumJ; nj++) {
                        arma::colvec ijInfo = iAsNN.col(nj);
                        int j = ijInfo(0);
                        int ksi = ijInfo(1);
                        double Bsjksi = Bs(j, ksi);
                        Vsi += std::pow(Bsjksi, 2) * (1 / Fs(j)) * KappaInv;
                    }
                }
            }
            Vsi = CholInv(Vsi);
            arma::mat cholVsi;
            try {
                cholVsi = arma::chol(Vsi);
            }
            catch (...) {
                cholVsi = getCholRobust(Vsi);
            }
            Vs.slice(i) = Vsi;
            cholVs.slice(i) = cholVsi;
        }
        int nnMaxNum;
        int nnMaxNumJ;
        for (arma::uword jj = 0; jj < K; jj++) {
            arma::mat AlphaJ = Alpha.slice(jj);
            arma::mat ZJ = Z.slice(jj);
            for (arma::uword lj = 0; lj < L - 1; lj++) {
                arma::mat ZJLMat = arma::reshape(ZJ.row(lj), O, M);
                arma::mat AlphaJLMat = arma::reshape(AlphaJ.row(lj), O, M);
                for (int i = 0; i < M; i++) {
                    arma::colvec muJLsi = ZJLMat.col(i);
                    if (i > 0) {
                        if (i < h) nnMaxNum = i - 1;
                        else nnMaxNum = h - 1;
                        arma::mat Bsi = Bs(arma::span(i, i), arma::span(0, nnMaxNum));
                        arma::uvec nni = nnInd(arma::span(i, i), arma::span(0, nnMaxNum)).t();
                        arma::mat AlphaJLNsi = arma::reshape(AlphaJLMat.cols(nni), (O * (nnMaxNum + 1)), 1);
                        muJLsi += arma::kron((1 / Fs(i)) * Bsi, KappaInv) * AlphaJLNsi;
                    }
                    if (i < M - 1) {
                        arma::mat iAsNN = whichJs(i, 0);
                        int iNumJ = size(iAsNN, 1);
                        if (!(iNumJ == 1 && iAsNN(0, 0) == 0)) {
                            for (arma::uword nj = 0; nj < iNumJ; nj++) {
                                arma::colvec ijInfo = iAsNN.col(nj);
                                int j = ijInfo(0);
                                if (j < h) nnMaxNumJ = j - 1;
                                else nnMaxNumJ = h - 1;
                                int ksi = ijInfo(1);
                                arma::colvec AlphaJLsj = AlphaJLMat.col(j);
                                arma::rowvec Bsj = Bs(arma::span(j, j), arma::span(0, nnMaxNumJ));
                                double Bsjksi = Bsj(ksi);
                                arma::uvec nnj = nnInd(arma::span(j, j), arma::span(0, nnMaxNumJ)).t();
                                arma::mat AlphaJLNsjMat = AlphaJLMat.cols(nnj);
                                AlphaJLNsjMat.col(ksi) = arma::colvec(O, arma::fill::zeros);
                                arma::mat AlphaJLNsj = arma::reshape(AlphaJLNsjMat, (O * (nnMaxNumJ + 1)), 1);
                                Bsj(ksi) = 0;
                                muJLsi += (Bsjksi * (1 / Fs(j)) * KappaInv) * (AlphaJLsj - arma::kron(Bsj, EyeO) * AlphaJLNsj);
                            }
                        }
                    }
                    arma::mat Vsi = Vs.slice(i);
                    arma::mat cholVsi = cholVs.slice(i);
                    AlphaJLMat.col(i) = rmvnormRcppNew(1, Vsi * muJLsi, cholVsi);
                }
                AlphaJ.row(lj) = arma::reshape(AlphaJLMat, 1, M * O);
            }
            Alpha.slice(jj) = AlphaJ;
        }
    }
          
    arma::cube logWeights = GetLogWeightsFixedL(Alpha, K, M, L, O);
    ParaCL.Alpha = Alpha;
    ParaCL.logWeights = logWeights;
    if (alphaWsFiles == 1 && storePostPara == true) {
        if (DatObj.spatPred == 1) {
            //std::string filename = "fixedLalphaIter" + to_string(keptIter) + ".txt";
            Alpha.save("fixedLalphaKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
        if (DatObj.storeW == 1) {
            logWeights.save("fixedLlogweightsKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
    }
    if (DatObj.storeW == 1 && storePostPara == true) {
        arma::cube Weights = GetWeightsFixedL(Alpha, K, M, L, O);
        if (alphaWsFiles == 1) {
            Weights.save("fixedLweightsKeptIter" + std::__cxx11::to_string(keptIter) + ".txt", arma::arma_ascii);
        }
        else { // if alphaWsFiles == 0
            ParaCL.Weights = Weights;
        }
    }
    return ParaCL;

}



//Function to sample zjl(s_i)'s using a Gibbs sampler step---------------------------------------------------------------
paraCLfixedL SampleZ(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL) {
  
  //Set data objects
  int K = DatObj.K;
  int M = DatObj.M;
  int O = DatObj.O;
  int L = DatObj.L;
  
  //Set parameter objects
  arma::cube Alpha = ParaCL.Alpha;
  arma::umat Xi = ParaCL.Xi;
  arma::cube Z = ParaCL.Z;
  
  //Declarations
  arma::uword Index1;
  arma::uword Index2;
  
  for (arma::uword j = 0; j < K; j++) {  
      arma::mat AlphaJ = Alpha.slice(j);
      arma::mat ZJ(L - 1, M * O);
      for (arma::uword i = 0; i < M; i++) {
          for (arma::uword o = 0; o < O; o++) {
              Index1 = i + M * o;
              Index2 = o + O * i;
              arma::uword XiOIJ = Xi(Index1, j);
              for (arma::uword l = 0; l < L - 1; l++) {
                  if (XiOIJ > l) ZJ(l, Index2) = rtnormRcppMSM(AlphaJ(l, Index2), 1, -arma::datum::inf, 0);
                  else if (XiOIJ == l) ZJ(l, Index2) = rtnormRcppMSM(AlphaJ(l, Index2), 1, 0, arma::datum::inf);
                  else ZJ(l, Index2) = arma::as_scalar(rnormRcpp(1, AlphaJ(l, Index2), 1));
              }
          }
      }
      Z.slice(j) = ZJ;   
  }
  
  //Update parameters object
  ParaCL.Z = Z;
  return ParaCL;
  
}



//Function to sample xi's using a Gibbs sampler step---------------------------------------------------------------
std::pair<VAR1paraFixedL, paraCLfixedL> SampleXi(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL) {
    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    int Nu = DatObj.Nu;
    arma::Col<int> SeqL = DatObj.SeqL;
    arma::mat YStarWide = DatObj.YStarWide;
    arma::mat EyeNu(Nu, Nu, arma::fill::eye);
    //Set parameter objects
    arma::umat Xi = ParaCL.Xi;
    arma::cube logWeights = ParaCL.logWeights;
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::cube Cov = Para.Cov;
    arma::mat Theta = ParaCL.Theta;
    arma::colvec Eta = Para.Eta;
    arma::mat Upsilon = Para.Upsilon;
    arma::colvec XBeta = Para.XBeta;
    arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);
    //Declarations
    arma::uword Index;
    for (arma::uword j = 0; j < K; j++) {
        arma::mat logWeightsJ = logWeights.slice(j);
        arma::colvec ThetaJ = Theta.col(j);
        for (arma::uword o = 0; o < O; o++) {
            arma::mat CovO = Cov(arma::span::all, arma::span(o, o), arma::span::all);
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                arma::colvec logWeightsOIJ = logWeightsJ.col(Index);
                arma::rowvec CovOI = CovO.row(i);
                arma::rowvec LambdaOI = Lambda.row(Index);
                arma::vec logProbsRaw(L);
                for (arma::uword l = 0; l < L; l++) {
                    //Obtain the un-normalized (raw) probabilities on the log scale
                    LambdaOI(j) = ThetaJ(l);
                    arma::rowvec Resid = YStarWide.row(Index) - XBetaMat.row(Index) - LambdaOI * BigPhi;
                    Resid = Resid % (1 / arma::sqrt(CovOI));
                    double ResidQ = arma::as_scalar(Resid * arma::trans(Resid));
                    double Likelihood = -0.5 * ResidQ;
                    logProbsRaw(l) = logWeightsOIJ(l) + Likelihood;
                }
                //Use log sum exponential trick to get normalized probabilities
                double Max = arma::max(logProbsRaw);
                double Delta = Max + log(arma::sum(arma::exp(logProbsRaw - Max)));
                arma::vec ProbsOIJ = arma::exp(logProbsRaw - Delta);
                //Sample a new label
                arma::uword XiOIJ = arma::as_scalar(sampleRcpp(SeqL, 1, true, ProbsOIJ)); 
                Xi(Index, j) = XiOIJ;
                //Update Lambda
                Lambda(Index, j) = Theta(XiOIJ, j);
            }
        }
    }
    //Update parameters object
    ParaCL.Xi = Xi;
    Para.Lambda = Lambda;
    return std::pair<VAR1paraFixedL, paraCLfixedL>(Para, ParaCL);
}



//Function to sample theta_jl's using a Gibbs sampler step---------------------------------------------------------------
std::pair<VAR1paraFixedL, paraCLfixedL> SampleTheta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL) {
    //Set data objects
    int K = DatObj.K;
    int L = DatObj.L;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat YStarWide = DatObj.YStarWide;
    arma::mat EyeNu(Nu, Nu, arma::fill::eye);
    //Set parameter objects
    arma::mat BigPhi = Para.BigPhi;
    arma::umat Xi = ParaCL.Xi;
    arma::mat Lambda = Para.Lambda;
    arma::cube Cov = Para.Cov;
    arma::mat CovMat = arma::reshape(arma::vectorise(Cov), M * O, Nu);
    arma::colvec Tau = ParaCL.Tau;
    arma::mat Theta = ParaCL.Theta;
    arma::colvec Eta = Para.Eta;
    arma::mat Upsilon = Para.Upsilon;
    arma::colvec XBeta = Para.XBeta;
    arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);
    for (arma::uword j = 0; j < K; j++) {
        arma::colvec EtaJ = arma::trans(BigPhi.row(j));
        arma::uvec XiJ = Xi.col(j);
        double TauJ = arma::as_scalar(Tau(j));
        arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
        LambdaMinusJ.shed_col(j);
        BigPhiMinusJ.shed_row(j);
        for (arma::uword l = 0; l < L; l++) {
            arma::uvec WhichJL = find(XiJ == l);
            int NJL = WhichJL.size();
            if (NJL == 0) { Theta(l, j) = arma::as_scalar(rnormRcpp(1, 0, sqrt(1 / TauJ))); }
            else {
                arma::mat CovInvL = 1 / CovMat.rows(WhichJL);
                arma::mat YJL = YStarWide.rows(WhichJL);
                arma::mat LambdaJL = LambdaMinusJ.rows(WhichJL);
                arma::mat XBetaComp = XBetaMat.rows(WhichJL);
                double ResidsJL = arma::as_scalar(arma::sum((YJL - XBetaComp - LambdaJL * BigPhiMinusJ) % CovInvL, 0) * EtaJ);
                double VarThetaJL = arma::as_scalar(1 / (arma::sum(CovInvL, 0) * (EtaJ % EtaJ) + TauJ));
                double MeanThetaJL = VarThetaJL * ResidsJL;
                double ThetaJL = arma::as_scalar(rnormRcpp(1, MeanThetaJL, sqrt(VarThetaJL)));
                Theta(l, j) = ThetaJL;
                arma::colvec LambdaJ = Lambda.col(j);
                arma::vec ThetaJLVec(1);
                ThetaJLVec(0) = ThetaJL;
                LambdaJ(WhichJL) = arma::repmat(ThetaJLVec, NJL, 1);
                Lambda.col(j) = LambdaJ;
            }            
        }
    }
    //Final updates
    arma::colvec Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
    //Update parameters object
    ParaCL.Theta = Theta;
    Para.Lambda = Lambda;
    Para.Mean = Mean;
    return std::pair<VAR1paraFixedL, paraCLfixedL>(Para, ParaCL);
}
