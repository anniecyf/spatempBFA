#define ARMA_DONT_PRINT_ERRORS //So the cholesky warning is suppressed
#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"

//Function to sample latent polya-gamma process using Gibbs sampling step------------------------------------------------
arma::mat SampleOmega(int f, VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para) {
  
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
VAR1datobjVaryLj SampleUpper(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, dataug DatAug) {
  
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
VAR1datobjVaryLj SampleLower(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, dataug DatAug) {
  
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
std::pair<VAR1datobjVaryLj, VAR1paraVaryLj> SampleY(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, dataug DatAug) {//called when if (any(FamilyInd != 0))
  
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
  return std::pair<VAR1datobjVaryLj, VAR1paraVaryLj>(DatObj, Para);

}



//Function to sample sigma2(s_i) using a Gibbs sampler step---------------------------------------------------------------
VAR1paraVaryLj SampleSigma2(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara) {
  
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
VAR1paraVaryLj SampleBeta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara) {
  
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
VAR1paraVaryLj SampleA(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara) {
  arma::mat A = Para.A;
  Para.A = A;
  return Para;
}



//Function to sample Upsilon using a Gibbs sampler step-------------------------------------------------------------------
VAR1paraVaryLj SampleUpsilon(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara) {
  
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
VAR1paraVaryLj SampleEta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara) {
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
  arma::colvec etaTminus1(K, arma::fill::zeros), etaTplus1(K), etaT(K), MeanEtaT(K);
  for (arma::uword t = 0; t < (Nu -  1); t++) { //for time = 1, 2, ..., (T - 1)
      if (t > 0) {etaTminus1 = BigPhi.col(t - 1);}
      etaTplus1 = BigPhi.col(t + 1);
      MeanEtaT = CovEta * (tLambdaSigmaInv * (YStarWide.col(t) - XBetaMat.col(t)) + UpsilonInv * A * etaTminus1 + arma::trans(A) * UpsilonInv * etaTplus1);
      etaT = rmvnormRcpp(1, MeanEtaT, CovEta);
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
VAR1paraVaryLj SampleDelta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara, listParaVaryLj ListPara) {  
    //Set data objects
    int K = DatObj.K;
    int GS = DatObj.GS;
    //Set parameter objects
    arma::field<arma::vec> Theta = ListPara.Theta;
    arma::colvec Delta = Para.Delta;
    arma::colvec LjVec = Para.LjVec;
    //Set hyperparameter objects
    double A1 = HyPara.A1;
    double A2 = HyPara.A2;
    //Gamma shrinkage prior
    if (GS == 1) {
        double AH;
        for (arma::uword h = 0; h < K; h++) {
            arma::colvec DeltaMinusH = Delta;
            DeltaMinusH(h) = 1;
            if (h == 0) AH = A1;
            else if (h > 0) AH = A2;
            double Resids = 0;
            for (arma::uword j = h; j < K; j++) {
                arma::colvec ThetaJ = Theta(j);
                double tThetaTheta = arma::as_scalar(arma::trans(ThetaJ) * ThetaJ);
                Resids += tThetaTheta * arma::prod(DeltaMinusH(arma::span(0, j)));
            }
            double Shape = AH;
            Shape += 0.5 * arma::sum(LjVec(arma::span(h, K - 1)));
            double Rate = 1 + 0.5 * Resids;
            Delta(h) = rgammaRcpp(Shape, Rate);
        }
    }
    //NO Gamma shrinkage prior
    else {
        for (arma::uword j = 0; j < K; j++) {
            int Lj = LjVec(j);
            arma::colvec ThetaJ = Theta(j);
            double tThetaTheta = arma::as_scalar(arma::trans(ThetaJ) * ThetaJ);          
            double Shape = A1 + 0.5 * Lj; 
            double Rate = A2 + 0.5 * tThetaTheta;
            Delta(j) = rgammaRcpp(Shape, Rate);
        }
    }
    //Update parameters object
    Para.Delta = Delta;
    if (GS == 1) { Para.Tau = arma::cumprod(Delta); }
    else { Para.Tau = Delta; } //if (GS == 0) 
    return Para; 
}



//Functions to sample new value of rho using a Metropolis sampler step-----------------------------------------------
std::pair<spatpara1VaryLj, VAR1metrobj> SampleRho(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj, listParaVaryLj ListPara) {
    int SpCorInd = DatObj.SpCorInd;
    if (SpCorInd == 0) {//spatial.structure = "continuous"
        //Set data objects
        int M = DatObj.M;
        int O = DatObj.O;
        arma::mat SpDist = DatObj.SpDist;
        int K = DatObj.K;
        arma::colvec ZeroOM(O * M, arma::fill::zeros);
        //Set parameter objects
        double Rho = SpatPara.Rho;
        arma::mat RootiKappa = Para.RootiKappa;
        arma::mat RootiSpCov = SpatPara.RootiSpCov;
        arma::colvec LjVec = Para.LjVec;
        arma::field<arma::mat> Alpha = ListPara.Alpha;

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
        RootiAlpha = arma::kron(RootiSpCov, RootiKappa);
        RootiAlphaProposal = arma::kron(RootiSpCovProposal, RootiKappa);
   
        //Loop over columns of j
        for (arma::uword j = 0; j < K; j++) {
            //Objects that depend on j
            int Lj = LjVec(j);
            if (Lj > 1) {
                arma::mat AlphaJ = Alpha(j);
                //Loop over mixture components
                for (arma::uword l = 0; l < Lj - 1; l++) {
                    arma::colvec AlphaJL = AlphaJ.row(l).t();
                    Component1A += lndMvn(AlphaJL, ZeroOM, RootiAlphaProposal);
                    Component1B += lndMvn(AlphaJL, ZeroOM, RootiAlpha);
                }
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
        int K = DatObj.K;
        arma::colvec ZeroOM(O * M, arma::fill::zeros); 
        //Set parameter objects
        double Rho = SpatPara.Rho;
        arma::mat KappaInv = Para.KappaInv;
        arma::mat RootiKappa = arma::trans(arma::chol(KappaInv)); //transpose of the cholesky factor of KappaInv instead; not the real rooti but corresponds to RootiSpCov here
        arma::mat RootiSpCov = SpatPara.RootiSpCov; //transpose of the cholesky factor of SpCovInv instead; not the real rooti but works here too
        arma::colvec LjVec = Para.LjVec;
        arma::field<arma::mat> Alpha = ListPara.Alpha;
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
        RootiAlpha = arma::kron(RootiSpCov, RootiKappa);
        RootiAlphaProposal = arma::kron(RootiSpCovProposal, RootiKappa);

        //Loop over columns of j
        for (arma::uword j = 0; j < K; j++) {
            //Objects that depend on j
            int Lj = LjVec(j);
            if (Lj > 1) {
                arma::mat AlphaJ = Alpha(j);
                //Loop over mixture components
                for (arma::uword l = 0; l < Lj - 1; l++) {
                    arma::colvec AlphaJL = AlphaJ.row(l).t();
                    Component1A += lndMvn(AlphaJL, ZeroOM, RootiAlphaProposal);
                    Component1B += lndMvn(AlphaJL, ZeroOM, RootiAlpha);
                }
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
    return std::pair<spatpara1VaryLj, VAR1metrobj>(SpatPara, MetrObj);
}

std::pair<spatpara2VaryLj, VAR1metrobj> SampleRho(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj, listParaVaryLj ListPara) {//when IS == 1 && SpCorInd == 0 && DatObj.SA == 1
    //Set data objects
    int M = DatObj.M;
    int O = DatObj.O;
    arma::mat SpDist = DatObj.SpDist;
    int K = DatObj.K;
    int h = DatObj.h;
    //Set parameter objects
    double Rho = SpatPara.Rho;
    arma::mat KappaInv = Para.KappaInv;
    arma::mat SpCovInvApprox = SpatPara.SpCovInv;
    double detSpCovApprox = SpatPara.detSpCov;
    arma::colvec LjVec = Para.LjVec;
    arma::field<arma::mat> Alpha = ListPara.Alpha;
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
    arma::mat kronInvDiff = arma::kron(SpCovInvDiff, KappaInv);
    arma::mat kronInvDiff2 = arma::kron(KappaInv, SpCovInvDiff);

    //Loop over columns of j
    for (arma::uword j = 0; j < K; j++) {
        //Objects that depend on j
        int Lj = LjVec(j);
        if (Lj > 1) {
            component1 -= logDetDiff * (Lj - 1);
            arma::mat AlphaJ = Alpha(j);
            //Loop over mixture components
            for (arma::uword l = 0; l < Lj - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                component1 -= arma::as_scalar(AlphaJL * kronInvDiff * arma::trans(AlphaJL));
            }
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
    }

    //Return output object
    return std::pair<spatpara2VaryLj, VAR1metrobj>(SpatPara, MetrObj);
}

std::pair<spatpara3, VAR1metrobj> SampleRho(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj, listParaVaryLj ListPara) {//when IS == 1 && SpCorInd == 0 && DatObj.SA == 1
    //Set data objects
    int M = DatObj.M;
    int O = DatObj.O;
    arma::mat SpDist = DatObj.SpDist;
    int K = DatObj.K;
    int h = DatObj.h;
    //Set parameter objects
    double Rho = SpatPara.Rho;
    arma::mat KappaInv = Para.KappaInv;
    arma::mat SpCovInvApprox = SpatPara.SpCovInv;
    double detSpCovApprox = SpatPara.detSpCov;
    arma::colvec LjVec = Para.LjVec;
    arma::field<arma::mat> Alpha = ListPara.Alpha;
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
    arma::mat kronInvDiff = arma::kron((SpCovInvApproxProposal - SpCovInvApprox), KappaInv);

    //Loop over columns of j
    for (arma::uword j = 0; j < K; j++) {
        //Objects that depend on j
        int Lj = LjVec(j);
        if (Lj > 1) {
            component1 -= logDetDiff * (Lj - 1);
            arma::mat AlphaJ = Alpha(j);
            //Loop over mixture components
            for (arma::uword l = 0; l < Lj - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                component1 -= arma::as_scalar(AlphaJL * kronInvDiff * arma::trans(AlphaJL));
            }
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
VAR1paraVaryLj SampleKappa(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, VAR1hypara HyPara, listParaVaryLj ListPara) {

    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;

    //Set parameter objects
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::colvec LjVec = Para.LjVec;
    arma::mat SpCovInv = SpatPara.SpCovInv;

    //Set hyperparameter objects
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;

    //Calculate moments
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj > 1) {
            arma::mat AlphaJ = Alpha(j);
            for (arma::uword l = 0; l < Lj - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                arma::mat AlphaJLMat = arma::reshape(AlphaJL, O, M);
                Resids += AlphaJLMat * SpCovInv * arma::trans(AlphaJLMat);
            }
        }
    }
    double Shape = SmallUpsilon + M * (arma::as_scalar(arma::sum(LjVec)) - K);
    arma::mat Rate = Resids + BigTheta;

    //Sample kappa
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

    //Update parameters object
    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;
  
}

VAR1paraVaryLj SampleKappa(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, VAR1hypara HyPara, listParaVaryLj ListPara) {

    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;

    //Set parameter objects
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::colvec LjVec = Para.LjVec;
    arma::mat SpCovInv = SpatPara.SpCovInv;

    //Set hyperparameter objects
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;

    //Calculate moments
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj > 1) {
            arma::mat AlphaJ = Alpha(j);
            for (arma::uword l = 0; l < Lj - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                arma::mat AlphaJLMat = arma::reshape(AlphaJL, O, M);
                Resids += AlphaJLMat * SpCovInv * arma::trans(AlphaJLMat);
            }
        }       
    }
    double Shape = SmallUpsilon + M * (arma::as_scalar(arma::sum(LjVec)) - K);
    arma::mat Rate = Resids + BigTheta;

    //Sample kappa
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

    //Update parameters object
    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;

}

VAR1paraVaryLj SampleKappa(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, VAR1hypara HyPara, listParaVaryLj ListPara) {

    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;

    //Set parameter objects
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::colvec LjVec = Para.LjVec;
    arma::mat SpCovInv = SpatPara.SpCovInv;

    //Set hyperparameter objects
    double SmallUpsilon = HyPara.SmallUpsilon;
    arma::mat BigTheta = HyPara.BigTheta;

    //Calculate moments
    arma::mat Resids(O, O, arma::fill::zeros);
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj > 1) {
            arma::mat AlphaJ = Alpha(j);
            for (arma::uword l = 0; l < Lj - 1; l++) {
                arma::rowvec AlphaJL = AlphaJ.row(l);
                arma::mat AlphaJLMat = arma::reshape(AlphaJL, O, M);
                Resids += AlphaJLMat * SpCovInv * arma::trans(AlphaJLMat);
            }
        }
    }
    double Shape = SmallUpsilon + M * (arma::as_scalar(arma::sum(LjVec)) - K);
    arma::mat Rate = Resids + BigTheta;

    //Sample kappa
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

    //Update parameters object
    Para.Kappa = Kappa;
    Para.KappaInv = KappaInv;
    Para.RootiKappa = RootiKappa;
    return Para;

}



//Functions to sample alphajl(s_i)'s using a Gibbs sampler step---------------------------------------------------------------
listParaVaryLj SampleAlpha(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    arma::colvec LjVec = Para.LjVec;
    arma::mat Kappa = Para.Kappa;   
    arma::mat SpCov = SpatPara.SpCov;
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::field<arma::mat> Weights = ListPara.Weights;
    arma::umat Xi = Para.Xi; // of dim MO X K
    arma::mat U = Para.U; // of dim MO X K
    arma::vec ZeroOM(O * M, arma::fill::zeros);
    arma::mat EyeOM(O * M, O * M, arma::fill::zeros);
    arma::mat CovAlpha = arma::kron(SpCov, Kappa);   
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj > 1) {
            arma::mat AlphaJ = Alpha(j); // of dim (Lj - 1) x MO, where Lj is the current MCMC iteration's Lj
            arma::mat WeightsJ = Weights(j); // of dim Lj x MO, where Lj is the current MCMC iteration's Lj
            arma::vec WjXiJOIsioVec(M * O);
            arma::uword Index1;
            arma::uword Index2;
            // prepare initial values (before the lj = 0 iteration)
            for (arma::uword o = 0; o < O; o++) {
                for (arma::uword i = 0; i < M; i++) {
                    Index1 = i + M * o;
                    arma::uword XiJOI = Xi(Index1, j);
                    WjXiJOIsioVec(Index1) = WeightsJ(XiJOI, Index1);
                }
            }
            // main iteration for the jth factor
            for (arma::uword lj = 0; lj < Lj - 1; lj++) {
                arma::rowvec AlphaJlj = AlphaJ.row(lj); // of length MO
                arma::vec upperVec(M * O); //upper truncation points for alphajlj's truncated multivariate normal distribution
                upperVec.fill(arma::datum::inf);
                arma::vec lowerVec(M * O); //lower truncation points for alphajlj's truncated multivariate normal distribution
                lowerVec.fill(-arma::datum::inf);
                for (arma::uword i = 0; i < M; i++) {
                    for (arma::uword o = 0; o < O; o++) {
                        Index1 = i + M * o;
                        Index2 = o + O * i;
                        double UJOI = U(Index1, j);
                        arma::uword XiJOI = Xi(Index1, j);
                        // if lj > XiJOI, then we don't need to change the default upper and lower values
                        if (lj == XiJOI) {
                            double denominator = WjXiJOIsioVec(Index1) / pnormRcpp(AlphaJlj(Index2));
                            lowerVec(Index2) = qnormRcpp(UJOI / denominator);
                        }
                        else if (lj < XiJOI) {
                            double denominator = WjXiJOIsioVec(Index1) / UpperpnormRcpp(AlphaJlj(Index2));
                            //Rcpp::Rcout << lj;//for debugging
                            upperVec(Index2) = UpperqnormRcpp(UJOI / denominator);
                        }
                    }
                }
                arma::rowvec AlphaJljNew = rtmvnormRcpp(ZeroOM, CovAlpha, lowerVec, upperVec);
                AlphaJ.row(lj) = AlphaJljNew;               
            }
            Alpha(j) = AlphaJ;
        }      
    }
    ListPara.Alpha = Alpha;
    ListPara.Weights = GetWeightsVaryLj(Alpha, K, M, O, LjVec);
    return ListPara;  
}

listParaVaryLj SampleAlpha(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    arma::colvec LjVec = Para.LjVec;
    arma::mat Kappa = Para.Kappa;
    arma::mat SpCov = SpatPara.SpCov;
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::field<arma::mat> Weights = ListPara.Weights;
    arma::umat Xi = Para.Xi; // of dim MO X K
    arma::mat U = Para.U; // of dim MO X K
    arma::vec ZeroOM(O * M, arma::fill::zeros);
    arma::mat EyeOM(O * M, O * M, arma::fill::zeros);
    arma::mat CovAlpha = arma::kron(SpCov, Kappa);
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj > 1) {
            arma::mat AlphaJ = Alpha(j); // of dim (Lj - 1) x MO, where Lj is the current MCMC iteration's Lj
            arma::mat WeightsJ = Weights(j); // of dim Lj x MO, where Lj is the current MCMC iteration's Lj
            arma::vec WjXiJOIsioVec(M * O);
            arma::uword Index1;
            arma::uword Index2;
            // prepare initial values (before the lj = 0 iteration)
            for (arma::uword o = 0; o < O; o++) {
                for (arma::uword i = 0; i < M; i++) {
                    Index1 = i + M * o;
                    arma::uword XiJOI = Xi(Index1, j);
                    WjXiJOIsioVec(Index1) = WeightsJ(XiJOI, Index1);
                }
            }
            // main iteration for the jth factor
            for (arma::uword lj = 0; lj < Lj - 1; lj++) {
                arma::rowvec AlphaJlj = AlphaJ.row(lj); // of length MO
                arma::vec upperVec(M * O); //upper truncation points for alphajlj's truncated multivariate normal distribution
                upperVec.fill(arma::datum::inf);
                arma::vec lowerVec(M * O); //lower truncation points for alphajlj's truncated multivariate normal distribution
                lowerVec.fill(-arma::datum::inf);
                for (arma::uword i = 0; i < M; i++) {
                    for (arma::uword o = 0; o < O; o++) {
                        Index1 = i + M * o;
                        Index2 = o + O * i;
                        double UJOI = U(Index1, j);
                        arma::uword XiJOI = Xi(Index1, j);
                        // if lj > XiJOI, then we don't need to change the default upper and lower values
                        if (lj == XiJOI) {
                            double denominator = WjXiJOIsioVec(Index1) / pnormRcpp(AlphaJlj(Index2));
                            lowerVec(Index2) = qnormRcpp(UJOI / denominator);
                        }
                        else if (lj < XiJOI) {
                            double denominator = WjXiJOIsioVec(Index1) / UpperpnormRcpp(AlphaJlj(Index2));
                            upperVec(Index2) = UpperqnormRcpp(UJOI / denominator);
                        }
                    }
                }
                arma::rowvec AlphaJljNew = rtmvnormRcpp(ZeroOM, CovAlpha, lowerVec, upperVec);
                AlphaJ.row(lj) = AlphaJljNew;
            }
            Alpha(j) = AlphaJ;
        }
    }
    ListPara.Alpha = Alpha;
    ListPara.Weights = GetWeightsVaryLj(Alpha, K, M, O, LjVec);
    return ListPara;
}

listParaVaryLj SampleAlpha(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int h = DatObj.h;
    arma::colvec LjVec = Para.LjVec;   
    arma::mat Bs = SpatPara.Bs;
    arma::colvec Fs = SpatPara.Fs;
    arma::umat nnInd = SpatPara.nnInd;
    arma::field<arma::mat> whichJs = SpatPara.whichJs;
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::field<arma::mat> Weights = ListPara.Weights;  
    arma::umat Xi = Para.Xi;// of dim MO X K
    arma::mat U = Para.U;// of dim MO X K

    if (O == 1) {
        arma::mat KappaInvMat = Para.KappaInv;
        double KappaInv = KappaInvMat(0, 0);
        arma::vec Vs(M);
        for (int i = 0; i < M; i++) {
            double Vsi = KappaInv / Fs(i);
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
            int Lj = LjVec(jj);
            if (Lj > 1) {
                arma::mat AlphaJ = Alpha(jj); // of dim (Lj - 1) x M, where Lj is the current MCMC iteration's Lj
                arma::mat WeightsJ = Weights(jj); // of dim Lj x M, where Lj is the current MCMC iteration's Lj 
                arma::vec WjXiJIsiVec(M);
                // prepare initial values (before the lj = 0 iteration)
                for (arma::uword i = 0; i < M; i++) {                    
                    arma::uword XiJI = Xi(i, jj);
                    WjXiJIsiVec(i) = WeightsJ(XiJI, i);
                }
                // main iteration for the jth factor                   
                for (arma::uword lj = 0; lj < Lj - 1; lj++) {
                    arma::colvec AlphaJlj = arma::trans(AlphaJ.row(lj));// of length M
                    for (int i = 0; i < M; i++) {
                        double muJljsi = 0;
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
                        arma::uword XiJI = Xi(i, jj);
                        double AlphaJljsiNew;
                        if (lj > XiJI) {
                            AlphaJljsiNew = arma::as_scalar(rnormRcpp(1, Vsi * muJljsi, sqrt(Vsi)));
                        }
                        else if (lj == XiJI) {
                            double denominator = WjXiJIsiVec(i) / pnormRcpp(AlphaJlj(i));
                            double UJI = U(i, jj);
                            double lower = qnormRcpp(UJI / denominator);
                            double upper = arma::datum::inf; //upper truncation point for alphajlj(si)'s truncated multivariate normal distribution
                            AlphaJljsiNew = rtnormRcppMSM(Vsi * muJljsi, sqrt(Vsi), lower, upper);
                        }
                        else { // if lj < XiJI
                            double denominator = WjXiJIsiVec(i) / UpperpnormRcpp(AlphaJlj(i));
                            double UJI = U(i, jj);
                            double upper = UpperqnormRcpp(UJI / denominator);
                            double lower = -arma::datum::inf; //lower truncation point for alphajlj(si)'s truncated multivariate normal distribution
                            AlphaJljsiNew = rtnormRcppMSM(Vsi * muJljsi, sqrt(Vsi), lower, upper);
                        }
                        AlphaJlj(i) = AlphaJljsiNew;
                    }
                    AlphaJ.row(lj) = arma::trans(AlphaJlj);
                }
                Alpha(jj) = AlphaJ;
            }
        }
    }
    else {//if O > 1
        arma::mat KappaInv = Para.KappaInv;
        arma::cube Vs(O, O, M);
        arma::mat EyeO(O, O, arma::fill::eye);
        for (int i = 0; i < M; i++) {
            arma::mat Vsi(O, O, arma::fill::zeros);
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
            Vs.slice(i) = Vsi;
        }
        int nnMaxNum;
        int nnMaxNumJ;
        for (arma::uword jj = 0; jj < K; jj++) {
            int Lj = LjVec(jj);
            if (Lj > 1) {
                arma::mat AlphaJ = Alpha(jj); // of dim (Lj - 1) x MO, where Lj is the current MCMC iteration's Lj
                arma::mat WeightsJ = Weights(jj); // of dim Lj x MO, where Lj is the current MCMC iteration's Lj  
                arma::vec WjXiJOIsioVec(M * O);
                arma::uword Index1;
                arma::uword Index2;
                // prepare initial values (before the lj = 0 iteration)
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index1 = i + M * o;
                        Index2 = o + O * i;
                        arma::uword XiJOI = Xi(Index1, jj);
                        WjXiJOIsioVec(Index2) = WeightsJ(XiJOI, Index1);// We need to use Index2 here for WjXiJOIsioVec instead for the sequential method to update Alpha
                    }
                }
                // main iteration for the jth factor
                for (arma::uword lj = 0; lj < Lj - 1; lj++) {
                    arma::mat AlphaJLMat = arma::reshape(AlphaJ.row(lj), O, M);
                    for (int i = 0; i < M; i++) {
                        arma::colvec muJLsi(O, arma::fill::zeros);
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
                        arma::colvec AlphaJLsi = AlphaJLMat.col(i); // of length O
                        arma::vec upperVec(O); //upper truncation points for alphajlj(si)'s truncated multivariate normal distribution
                        upperVec.fill(arma::datum::inf);
                        arma::vec lowerVec(O); //lower truncation points for alphajlj(si)'s truncated multivariate normal distribution
                        lowerVec.fill(-arma::datum::inf);
                        arma::vec WjXiJIsiVec = WjXiJOIsioVec(arma::span(i * O, (i + 1) * O - 1)); // of length O
                        bool truncated = false;
                        for (arma::uword o = 0; o < O; o++) {
                            Index1 = i + M * o;
                            double UJOI = U(Index1, jj);
                            arma::uword XiJOI = Xi(Index1, jj);
                            // if lj > XiJOI, then we don't need to change the default upper and lower values
                            if (lj == XiJOI) {
                                double denominator = WjXiJIsiVec(o) / pnormRcpp(AlphaJLsi(o));
                                lowerVec(o) = qnormRcpp(UJOI / denominator);
                                truncated = true;
                            }
                            else if (lj < XiJOI) {
                                double denominator = WjXiJIsiVec(o) / UpperpnormRcpp(AlphaJLsi(o));
                                upperVec(o) = UpperqnormRcpp(UJOI / denominator);
                                truncated = true;
                            }                            
                        }
                        arma::colvec AlphaJLsiNew;
                        if (truncated) { AlphaJLsiNew = arma::trans(rtmvnormRcpp(Vsi * muJLsi, Vsi, lowerVec, upperVec)); }
                        else { AlphaJLsiNew = rmvnormRcpp(1, Vsi * muJLsi, Vsi); }
                        AlphaJLMat.col(i) = AlphaJLsiNew;                       
                    }
                    AlphaJ.row(lj) = arma::reshape(AlphaJLMat, 1, M * O);
                }
                Alpha(jj) = AlphaJ;
            }            
        }
    }
    ListPara.Alpha = Alpha;
    ListPara.Weights = GetWeightsVaryLj(Alpha, K, M, O, LjVec);
    return ListPara;
}



//Function to sample xi's using a Gibbs sampler step---------------------------------------------------------------
std::pair<VAR1paraVaryLj, listParaVaryLj> SampleXi(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat YStarWide = DatObj.YStarWide;
    arma::mat EyeNu(Nu, Nu, arma::fill::eye);
    arma::umat Xi = Para.Xi;
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::cube Cov = Para.Cov;
    arma::field<arma::vec> Theta = ListPara.Theta;   
    arma::colvec Eta = Para.Eta;
    arma::mat Upsilon = Para.Upsilon;
    arma::colvec LjVecOld = Para.LjVec;
    arma::mat U = Para.U;
    arma::colvec XBeta = Para.XBeta;
    arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);
    arma::colvec LjVec = UpdateLjVec(U, ListPara.Weights, K, M, O, LjVecOld);
    arma::field<arma::mat> Alpha = UpdateAlpha(ListPara.Alpha, K, M, O, LjVec);//for each j, the new AlphaJ is a submatrix of AlphaJOld (its first Lj.new - 1 rows)
    arma::field<arma::mat> Weights = GetWeightsVaryLj(Alpha, K, M, O, LjVec);
    arma::uword Index;
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        arma::vec ThetaJ = Theta(j);
        if (Lj == 1) {
            arma::uvec XiJ(M * O, arma::fill::zeros);
            Xi.col(j) = XiJ;
            arma::colvec LambdaJ(M * O, arma::fill::zeros);
            LambdaJ += ThetaJ(0);
            Lambda.col(j) = LambdaJ;
        }
        else {
            arma::mat WeightsJ = Weights(j);
            //WeightsJ = WeightsJ(arma::span(0, Lj - 1), arma::span::all);
            arma::colvec UJ = U.col(j);
            for (arma::uword o = 0; o < O; o++) {
                arma::mat CovO = Cov(arma::span::all, arma::span(o, o), arma::span::all);
                for (arma::uword i = 0; i < M; i++) {
                    Index = i + M * o;
                    arma::colvec WeightsOIJ = WeightsJ.col(Index);
                    arma::rowvec CovOI = CovO.row(i);
                    arma::rowvec LambdaOI = Lambda.row(Index);
                    double UOIJ = UJ(Index);
                    //arma::vec ProbsOIJ(Lj, arma::fill::zeros);
                    arma::uvec whichPosProb = arma::find(WeightsOIJ > UOIJ);//at least one -- ensured by our algorithm
                    int posProbSize = whichPosProb.size();
                    arma::vec logProbsRaw(posProbSize);
                    //Loop over clusters l
                    for (int linclu = 0; linclu < posProbSize; linclu++) {
                        arma::uword l = whichPosProb(linclu);
                        LambdaOI(j) = ThetaJ(l);
                        arma::rowvec Resid = YStarWide.row(Index) - XBetaMat.row(Index) - LambdaOI * BigPhi;
                        Resid = Resid % (1 / arma::sqrt(CovOI));
                        double ResidQ = arma::as_scalar(Resid * arma::trans(Resid));
                        double logLikelihood = -0.5 * ResidQ;
                        logProbsRaw(linclu) = logLikelihood;
                        //ProbsOIJ(l) = exp(logLikelihood);
                    }
                    //Use log sum exponential trick to get normalized probabilities
                    double Max = arma::max(logProbsRaw);
                    double Delta = Max + log(arma::sum(arma::exp(logProbsRaw - Max)));
                    arma::vec ProbsOIJpos = arma::exp(logProbsRaw - Delta);
                    //ProbsOIJ(whichPosProb) = ProbsOIJpos;
                    //double sumProbProp = arma::sum(ProbsOIJ);
                    //ProbsOIJ /= sumProbProp;               
                    //arma::Col<int> SeqLjVec(Lj);
                    //for (arma::uword seq = 0; seq < Lj; seq++) SeqLjVec(seq) = seq;
                    //XiOIJ = arma::as_scalar(sampleRcpp(SeqLjVec, 1, true, ProbsOIJ));
                    arma::Col<int> SeqLjVec(posProbSize);
                    for (arma::uword seq = 0; seq < posProbSize; seq++) SeqLjVec(seq) = whichPosProb(seq);
                    arma::uword XiOIJ;
                    XiOIJ = arma::as_scalar(sampleRcpp(SeqLjVec, 1, true, ProbsOIJpos));
                    Xi(Index, j) = XiOIJ;
                    Lambda(Index, j) = ThetaJ(XiOIJ);
                }
            }
        }       
    }
    Para.LjVec = LjVec;
    ListPara.Weights = Weights;
    ListPara.Alpha = Alpha;
    Para.Xi = Xi;
    Para.Lambda = Lambda;
    return std::pair<VAR1paraVaryLj, listParaVaryLj>(Para, ListPara);
}



//Function to sample theta_jl's using a Gibbs sampler step---------------------------------------------------------------
std::pair<VAR1paraVaryLj, listParaVaryLj> SampleTheta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat YStarWide = DatObj.YStarWide;
    arma::mat EyeNu(Nu, Nu, arma::fill::eye);
    arma::mat BigPhi = Para.BigPhi;
    arma::umat Xi = Para.Xi;
    arma::mat Lambda = Para.Lambda;
    arma::cube Cov = Para.Cov;
    arma::mat CovMat = arma::reshape(arma::vectorise(Cov), M * O, Nu);
    arma::vec Tau = Para.Tau;
    arma::field<arma::vec> Theta(K);
    arma::colvec Eta = Para.Eta;
    arma::colvec LjVec = Para.LjVec;
    arma::mat Upsilon = Para.Upsilon;
    arma::colvec XBeta = Para.XBeta;
    arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        arma::vec ThetaJ(Lj);
        arma::colvec EtaJ = arma::trans(BigPhi.row(j));
        arma::uvec XiJ = Xi.col(j);
        double TauJ = arma::as_scalar(Tau.row(j));
        arma::mat LambdaMinusJ = Lambda, BigPhiMinusJ = BigPhi;
        LambdaMinusJ.shed_col(j);
        BigPhiMinusJ.shed_row(j);  
        for (arma::uword l = 0; l < Lj; l++) {
            arma::uvec WhichJL = find(XiJ == l);
            int NJL = WhichJL.size();
            double ThetaJL;
            if (NJL == 0) { ThetaJ(l) = arma::as_scalar(rnormRcpp(1, 0, sqrt(1 / TauJ))); }
            else {
                arma::mat CovInvL = 1 / CovMat.rows(WhichJL);
                arma::mat YJL = YStarWide.rows(WhichJL);
                arma::mat LambdaJL = LambdaMinusJ.rows(WhichJL);
                arma::mat XBetaComp = XBetaMat.rows(WhichJL);
                double ResidsJL = arma::as_scalar(arma::sum((YJL - XBetaComp - LambdaJL * BigPhiMinusJ) % CovInvL, 0) * EtaJ);
                double VarThetaJL = arma::as_scalar(1 / (arma::sum(CovInvL, 0) * (EtaJ % EtaJ) + TauJ));
                double MeanThetaJL = VarThetaJL * ResidsJL;
                ThetaJL = arma::as_scalar(rnormRcpp(1, MeanThetaJL, sqrt(VarThetaJL)));
                ThetaJ(l) = ThetaJL;
                arma::colvec LambdaJ = Lambda.col(j);
                arma::vec ThetaJLVec(1);
                ThetaJLVec(0) = ThetaJL;
                LambdaJ(WhichJL) = arma::repmat(ThetaJLVec, NJL, 1);
                Lambda.col(j) = LambdaJ;
            }
        }
        Theta(j) = ThetaJ;
    }

    //Final updates
    arma::colvec Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;

    //Update parameters object
    ListPara.Theta = Theta;
    Para.Lambda = Lambda;
    Para.Mean = Mean;
    return std::pair<VAR1paraVaryLj, listParaVaryLj>(Para, ListPara);
}



//Function to sample Latent U using a Gibbs sampler step-------------------------------------------------------------------
VAR1paraVaryLj SampleU(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, listParaVaryLj ListPara) {
    //Set data objects
    int M = DatObj.M;
    int K = DatObj.K;
    int O = DatObj.O;

    //Set parameters
    arma::field<arma::mat> Weights = ListPara.Weights;
    arma::umat Xi = Para.Xi;
    arma::colvec LjVec = Para.LjVec;

    //Sample latent U
    arma::mat U(M * O, K);
    arma::uword Index;
    for (arma::uword j = 0; j < K; j++) {
        arma::mat WeightsJ = Weights(j);
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                double WeightsOIJ = WeightsJ(Xi(Index, j), Index);
                U(Index, j) = randuRcpp() * WeightsOIJ;
            }
        }
    }
    Para.U = U;
    return Para;
}
