#define ARMA_DONT_PRINT_ERRORS //So the cholesky warning is suppressed
#include <RcppArmadillo.h>
#include "MCMC_bfaSpatemp.h"

//Function to sample latent polya-gamma process using Gibbs sampling step------------------------------------------------
arma::mat SampleOmega(int f, datobjVaryLj DatObj, paraVaryLj Para) {
  
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
datobjVaryLj SampleUpper(datobjVaryLj DatObj, paraVaryLj Para, dataug DatAug) {
  
  //Set data objects
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
datobjVaryLj SampleLower(datobjVaryLj DatObj, paraVaryLj Para, dataug DatAug) {
  
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
std::pair<datobjVaryLj, paraVaryLj> SampleY(datobjVaryLj DatObj, paraVaryLj Para, dataug DatAug) {//called when if (any(FamilyInd != 0))
  
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
  return std::pair<datobjVaryLj, paraVaryLj>(DatObj, Para);

}



//Function to sample sigma2(s_i) using a Gibbs sampler step---------------------------------------------------------------
paraVaryLj SampleSigma2(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara) {
  
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



//Function to sample new value of psi using a Metropolis sampler step-----------------------------------------------
std::pair<paraVaryLj, metrobj> SamplePsi(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara, metrobj MetrObj) {
  //Set data objects
  arma::mat TimeDist = DatObj.TimeDist;
  int Nu = DatObj.Nu;
  int K = DatObj.K;
  int TempCorInd = DatObj.TempCorInd;
  int seasonPeriod = DatObj.seasonPeriod;
  arma::colvec ZeroKbyNu(K * Nu, arma::fill::zeros); 
  int ET = DatObj.ET;
  //Set parameter objects
  double Psi = Para.Psi;
  arma::mat RootiHPsi = Para.RootiHPsi;
  arma::mat RootiUpsilon = Para.RootiUpsilon;
  arma::colvec Eta = Para.Eta;
  //Set hyperparameter objects
  double APsi = HyPara.APsi;
  double BPsi = HyPara.BPsi;
  double Gamma = HyPara.Gamma;
  double Beta = HyPara.Beta;
  //Set metropolis objects
  double MetropPsi = sqrt(MetrObj.MetropPsi);
  double AcceptancePsi = MetrObj.AcceptancePsi;
  //Transform current state to real line
  double BigDelta = log((Psi - APsi) / (BPsi - Psi));
  
  double PsiProposal, BigDeltaProposal;
  arma::mat RootiHPsiProposal(Nu, Nu);
  arma::mat CholHPsiProposal(Nu, Nu), HPsiProposal(Nu, Nu);//only used when ET==1
  //double psi = Psi;//default AR(1), SAR(1)
  //if (TempCorInd == 0 || TempCorInd == 3) {psi = exp(-Psi);}// for "exponential" or "sexponential" temporal.structure
  if (ET == 1) {// when equalTimeDist = TRUE
    //Sample a new Proposal
    BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropPsi));
    //Compute Phi Proposal
    PsiProposal = (BPsi * exp(BigDeltaProposal) + APsi) / (1 + exp(BigDeltaProposal));
    RootiHPsiProposal = getRootiH(Psi, Nu, TempCorInd, seasonPeriod);
  }
  else{//if ET == 0, i.e., when equalTimeDist = FALSE
    //Numerical fix for when the propopsal cholesky doesn't exist
    bool Cholesky = false;
    while (!Cholesky) {
      //Sample a new Proposal
      BigDeltaProposal = arma::as_scalar(rnormRcpp(1, BigDelta, MetropPsi));
      //Compute Phi Proposal
      PsiProposal = (BPsi * exp(BigDeltaProposal) + APsi) / (1 + exp(BigDeltaProposal));
      //Proposal temporal correlation
      HPsiProposal = getH(PsiProposal, TempCorInd, TimeDist, Nu, seasonPeriod);
      Cholesky = arma::chol(CholHPsiProposal, HPsiProposal); 
    }
    RootiHPsiProposal = arma::inv(arma::trimatu(CholHPsiProposal));
  }
  
  //Eta structure components
  arma::mat RootiEta = arma::kron(RootiHPsi, RootiUpsilon);
  arma::mat RootiEtaProposal = arma::kron(RootiHPsiProposal, RootiUpsilon);
  
  double Component1A = lndMvn(Eta, ZeroKbyNu, RootiEtaProposal);
  double Component1B = lndMvn(Eta, ZeroKbyNu, RootiEta);
  double Component1 = Component1A - Component1B; 
  //Prior component
  double Component2 = 0; //exponential
  if ((TempCorInd == 1) || ((TempCorInd == 2))) { //ar1 or sar1 transformed beta prior for psi
    double Component2A = (Gamma - 1) * log(1 + PsiProposal) + (Beta - 1) * log(1 - PsiProposal);
    double Component2B = (Gamma - 1) * log(1 + Psi) + (Beta - 1) * log(1 - Psi);
    Component2 = Component2A - Component2B;
  }
  double Component3A = BigDeltaProposal;
  double Component3B = BigDelta;
  double Component3 = Component3A - Component3B;
  double Component4 = 2 * log((1 + exp(BigDelta)) / (1 + exp(BigDeltaProposal)));
  //Compute log acceptance ratio
  double LogR = Component1 + Component2 + Component3 + Component4;
  
  //Metropolis update
  double RandU = randuRcpp();
  if (log(RandU) < LogR) {
    //Keep Count of Acceptances
    AcceptancePsi++;
    MetrObj.AcceptancePsi = AcceptancePsi;
    //Update parameters object
    Para.RootiHPsi = RootiHPsiProposal;
    Para.Psi = PsiProposal;
    if (ET == 1) {
      Para.HPsiInv = getInvH(Psi, Nu, TempCorInd, seasonPeriod);
    }
    else {//if ET == 0
      Para.HPsi = HPsiProposal;
      Para.HPsiInv = RootiHPsiProposal * arma::trans(RootiHPsiProposal);
    }
  }
  
  //Return output object
  return std::pair<paraVaryLj, metrobj>(Para, MetrObj);
}



//Function to sample Upsilon using a Gibbs sampler step-------------------------------------------------------------------
paraVaryLj SampleUpsilon(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara) {
  
  //Set data objects
  int Nu = DatObj.Nu;
  int K = DatObj.K;

  //Set parameters
  arma::mat BigPhi = Para.BigPhi;
  arma::mat HPsiInv = Para.HPsiInv;

  //Set hyperparameter objects
  double Zeta = HyPara.Zeta;
  arma::mat Omega = HyPara.Omega;
  
  //Compute SPhiPsi
  arma::mat SPhiPsi = BigPhi * HPsiInv * arma::trans(BigPhi);
  
  //Sample Upsilon
  double n = Zeta + Nu;
  arma::mat V = SPhiPsi + Omega;
  arma::mat Upsilon(K, K), UpsilonInv(K, K), RootiUpsilon(K, K);
  if (K > 1) {
    UpsilonInv = rwishRcpp(n, CholInv(V));
    Upsilon = CholInv(UpsilonInv);
    RootiUpsilon = GetRooti(Upsilon);
  } 
  else {
    Upsilon(0, 0) = rigammaRcpp(0.5 * n, 0.5 * arma::as_scalar(V));
    UpsilonInv = 1 / Upsilon;
    RootiUpsilon = arma::sqrt(UpsilonInv);
  }
  
  //Update parameters object
  Para.Upsilon = Upsilon;
  Para.UpsilonInv = UpsilonInv;
  Para.RootiUpsilon = RootiUpsilon;
  return Para;
}



//Function to sample beta using a Gibbs sampler step---------------------------------------------------------------
paraVaryLj SampleBeta(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara) {
  
  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  arma::mat YStarWide = DatObj.YStarWide;
  int Nu = DatObj.Nu;
  int P = DatObj.P;
  arma::mat X = DatObj.X;
  arma::Col<int> Indeces = DatObj.Indeces;
  
  //Set parameters
  arma::mat Lambda = Para.Lambda;
  //arma::mat Sigma2 = Para.Sigma2;
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



//Function to sample eta using a Gibbs sampler step---------------------------------------------------------------
paraVaryLj SampleEta(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara) {
  //Set data objects
  arma::mat EyeNu = DatObj.EyeNu;
  int K = DatObj.K;
  arma::mat EyeK(K, K, arma::fill::eye);
  arma::mat YStarWide = DatObj.YStarWide; 
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  int O = DatObj.O;
  int ET = DatObj.ET;
  int TempCorInd = DatObj.TempCorInd;
  int seasonPeriod = DatObj.seasonPeriod;
  int IT = DatObj.IT;
      
  //Set parameters
  arma::mat BigPhi = Para.BigPhi;
  arma::mat Lambda = Para.Lambda;
  arma::cube Cov = Para.Cov;
  arma::mat UpsilonInv = Para.UpsilonInv;
  arma::colvec XBeta = Para.XBeta;
  arma::mat XBetaMat = arma::reshape(XBeta, M * O, Nu);
  
  if (IT == 0) {
      arma::mat CondPrecEta = UpsilonInv;
      arma::colvec EtaT(K);
      arma::mat SigmaInv = arma::diagmat(arma::vectorise(1 / Cov.slice(0)));
      arma::mat tLambdaSigmaInv = arma::trans(Lambda) * SigmaInv;
      arma::mat CovEtaT = CholInv(tLambdaSigmaInv * Lambda + CondPrecEta);
      arma::mat cholSigma;
      try {
          cholSigma = arma::chol(CovEtaT);
      }
      catch (...) {
          cholSigma = getCholRobust(CovEtaT);
      }
      //Loop over t
      for (arma::uword t = 0; t < Nu; t++) {
          //Sample EtaT
          arma::colvec MeanEtaT = CovEtaT * (tLambdaSigmaInv * (YStarWide.col(t) - XBetaMat.col(t)));
          EtaT = rmvnormRcppNew(1, MeanEtaT, cholSigma);
          BigPhi.col(t) = EtaT;
          //End loop over t
      }
  }
  else {// if IT==1
      if (ET == 1) {//if equalTimeDist = TRUE
          double Psi = Para.Psi;
          double psi = Psi;//default ar1 or sar 1
          if (TempCorInd == 0 || TempCorInd == 3) {//exponential or sexponential
              psi = exp(-Psi);
          }
          double HStarInv = (1 + std::pow(psi, 2)) / (1 - std::pow(psi, 2));
          //for t=2,3,...,(T-1); HStarInv = 1 / (1 - psi ^ 2) for t=1 and t=T (if no seasonality);
          //for t=d+1,d+2,...,(T-d); HStarInv = 1 / (1 - psi ^ 2) for t=1,...,d, and t=T-d+1,...,T (if with seasonality period = d);
          //in C++ we should minus 1 to the actual t
          arma::mat CondPrecEta = HStarInv * UpsilonInv;
          arma::colvec EtaT(K);
          arma::mat SigmaInv = arma::diagmat(arma::vectorise(1 / Cov.slice(0)));
          arma::mat tLambdaSigmaInv = arma::trans(Lambda) * SigmaInv;
          //Loop over t
          for (arma::uword t = 0; t < Nu; t++) {
              arma::mat BigPhiMinusT = BigPhi;
              BigPhiMinusT.shed_col(t);
              arma::rowvec HPlus((Nu - 1), arma::fill::zeros);
              if (TempCorInd == 0 || TempCorInd == 1) {//no seasonality
                  if (t == 0) {
                      HPlus(0) = psi;
                      HStarInv = 1 / (1 - std::pow(psi, 2));
                      CondPrecEta = HStarInv * UpsilonInv;
                  }
                  else if (t == (Nu - 1)) {
                      HPlus((Nu - 2)) = psi;
                      HStarInv = 1 / (1 - std::pow(psi, 2));
                      CondPrecEta = HStarInv * UpsilonInv;
                  }
                  else {
                      HPlus(t) = psi / (1 + std::pow(psi, 2));
                      HPlus(t - 1) = HPlus(t);
                  }
              }
              else if (TempCorInd == 2 || TempCorInd == 3) {//with seasonality
                  if (t < seasonPeriod) {
                      HStarInv = 1 / (1 - std::pow(psi, 2));
                      CondPrecEta = HStarInv * UpsilonInv;
                      HPlus(seasonPeriod + t - 1) = psi;
                  }
                  else if (t >= (Nu - seasonPeriod)) {
                      HStarInv = 1 / (1 - std::pow(psi, 2));
                      CondPrecEta = HStarInv * UpsilonInv;
                      HPlus((t - seasonPeriod)) = psi;
                  }
                  else {
                      HPlus(t - seasonPeriod) = psi / (1 + std::pow(psi, 2));
                      HPlus(t + seasonPeriod - 1) = HPlus(t - seasonPeriod);
                  }
              }

              arma::colvec CondMuEta = arma::kron(HPlus, EyeK) * arma::vectorise(BigPhiMinusT);
              //Sample EtaT
              arma::mat CovEtaT = CholInv(tLambdaSigmaInv * Lambda + CondPrecEta);
              arma::colvec MeanEtaT = CovEtaT * (tLambdaSigmaInv * (YStarWide.col(t) - XBetaMat.col(t)) + CondPrecEta * CondMuEta);
              EtaT = rmvnormRcpp(1, MeanEtaT, CovEtaT);
              BigPhi.col(t) = EtaT;
              //End loop over t
          }
      }
      else if (ET == 0) {//if ET = 0, i.e., if equalTimeDist = FALSE
        //Declarations
          arma::mat HPsi = Para.HPsi;
          arma::vec SeqNu = arma::linspace<arma::vec>(0, Nu - 1, Nu);
          arma::uvec IndecesT(1), IndecesMinusT(Nu - 1);
          arma::colvec EtaT(K);
          arma::mat SigmaInv = arma::diagmat(arma::vectorise(1 / Cov.slice(0)));
          arma::mat tLambdaSigmaInv = arma::trans(Lambda) * SigmaInv;
          //Loop over t
          for (arma::uword t = 0; t < Nu; t++) {
              //Conditional moments
              IndecesMinusT = find(SeqNu != t);
              IndecesT(0) = t;
              arma::mat BigPhiMinusT = BigPhi;
              BigPhiMinusT.shed_col(t);
              arma::rowvec HPlus = HPsi(IndecesT, IndecesMinusT) * CholInv(HPsi(IndecesMinusT, IndecesMinusT));
              double HStarInv = 1 / (arma::as_scalar(HPsi(IndecesT, IndecesT) - HPlus * HPsi(IndecesMinusT, IndecesT)));
              arma::colvec CondMuEta = arma::kron(HPlus, EyeK) * arma::vectorise(BigPhiMinusT);
              arma::mat CondPrecEta = HStarInv * UpsilonInv;
              //Sample EtaT
              arma::mat CovEtaT = CholInv(tLambdaSigmaInv * Lambda + CondPrecEta);
              arma::colvec MeanEtaT = CovEtaT * (tLambdaSigmaInv * (YStarWide.col(t) - XBetaMat.col(t)) + CondPrecEta * CondMuEta);
              EtaT = rmvnormRcpp(1, MeanEtaT, CovEtaT);
              BigPhi.col(t) = EtaT;
              //End loop over visits
          }
      }
  }
  
  //Update parameters dependent on eta
  arma::colvec Eta = arma::vectorise(BigPhi);  
  //Update parameters object
  Para.Eta = Eta;
  Para.BigPhi = BigPhi;
  Para.Mean = arma::kron(EyeNu, Lambda) * Eta + XBeta;
  return Para;
}



//Function to sample delta using a Gibbs sampler step---------------------------------------------------------------
paraVaryLj SampleDelta(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara, listParaVaryLj ListPara) {  
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
std::pair<spatpara1VaryLj, metrobj> SampleRho(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, hypara HyPara, metrobj MetrObj, listParaVaryLj ListPara) {
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
    return std::pair<spatpara1VaryLj, metrobj>(SpatPara, MetrObj);
}

std::pair<spatpara2VaryLj, metrobj> SampleRho(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, hypara HyPara, metrobj MetrObj, listParaVaryLj ListPara) {//when IS == 1 && SpCorInd == 0 && DatObj.SA == 1
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
    arma::sp_mat kronInvDiff(arma::kron(SpCovInvDiff, KappaInv));
    //arma::mat kronInvDiff = arma::kron(SpCovInvDiff, KappaInv);

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
    return std::pair<spatpara2VaryLj, metrobj>(SpatPara, MetrObj);
}

std::pair<spatpara3, metrobj> SampleRho(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, hypara HyPara, metrobj MetrObj, listParaVaryLj ListPara) {//when IS == 1 && SpCorInd == 0 && DatObj.SA == 1
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
    arma::sp_mat kronInvDiff(arma::kron((SpCovInvApproxProposal - SpCovInvApprox), KappaInv));
    //arma::mat kronInvDiff = arma::kron((SpCovInvApproxProposal - SpCovInvApprox), KappaInv);

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
    return std::pair<spatpara3, metrobj>(SpatPara, MetrObj);
}





//Functions to sample kappa using a Gibbs sampler step---------------------------------------------------------------
paraVaryLj SampleKappa(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, hypara HyPara, listParaVaryLj ListPara) {

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

paraVaryLj SampleKappa(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, hypara HyPara, listParaVaryLj ListPara) {

    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;

    //Set parameter objects
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::colvec LjVec = Para.LjVec;
    arma::sp_mat SpCovInv(SpatPara.SpCovInv);
    //arma::mat SpCovInv = SpatPara.SpCovInv;

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

paraVaryLj SampleKappa(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, hypara HyPara, listParaVaryLj ListPara) {

    //Set data objects
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;

    //Set parameter objects
    arma::field<arma::mat> Alpha = ListPara.Alpha;
    arma::colvec LjVec = Para.LjVec;
    arma::sp_mat SpCovInv(SpatPara.SpCovInv);
    //arma::mat SpCovInv = SpatPara.SpCovInv;

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
listParaVaryLj SampleAlpha(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, listParaVaryLj ListPara) {
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

listParaVaryLj SampleAlpha(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, listParaVaryLj ListPara) {
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

listParaVaryLj SampleAlpha(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, listParaVaryLj ListPara) {
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
            double Vsi = 1 / Fs(i);
            if (i < M - 1) {
                arma::mat iAsNN = whichJs(i, 0);
                int iNumJ = size(iAsNN, 1);
                if (!(iNumJ == 1 && iAsNN(0, 0) == 0)) {
                    for (arma::uword nj = 0; nj < iNumJ; nj++) {
                        arma::colvec ijInfo = iAsNN.col(nj);
                        int j = ijInfo(0);
                        int ksi = ijInfo(1);
                        double Bsjksi = Bs(j, ksi);
                        Vsi += std::pow(Bsjksi, 2) / Fs(j);
                    }
                }
            }
            Vsi = 1 / Vsi;
            Vs(i) = Vsi / KappaInv;
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
                            muJljsi += arma::as_scalar(Bsi * AlphaJljNsi) / Fs(i);
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
                                    muJljsi += Bsjksi / Fs(j) * (AlphaJljsj - arma::as_scalar(Bsj * AlphaJljNsj));
                                }
                            }
                        }
                        muJljsi *= KappaInv;
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
        arma::mat EyeO(O, O, arma::fill::eye);
        arma::mat KappaInv = Para.KappaInv;
        arma::mat Kappa = Para.Kappa;
        arma::mat cholKappa;
        try {
            cholKappa = arma::chol(Kappa);
        }
        catch (...) {
            cholKappa = getCholRobust(Kappa);
        }
        arma::cube Vs(O, O, M);   
        arma::cube cholVs(O, O, M);
        for (int i = 0; i < M; i++) {
            double VsiSca = 1 / Fs(i);
            if (i < M - 1) {
                arma::mat iAsNN = whichJs(i, 0);
                int iNumJ = size(iAsNN, 1);
                if (!(iNumJ == 1 && iAsNN(0, 0) == 0)) {
                    for (arma::uword nj = 0; nj < iNumJ; nj++) {
                        arma::colvec ijInfo = iAsNN.col(nj);
                        int j = ijInfo(0);
                        int ksi = ijInfo(1);
                        double Bsjksi = Bs(j, ksi);
                        VsiSca += std::pow(Bsjksi, 2) / Fs(j);
                    }
                }
            }
            arma::mat Vsi = (1 / VsiSca) * Kappa;
            arma::mat cholVsi = (1 / sqrt(VsiSca)) * cholKappa;
            Vs.slice(i) = Vsi;
            cholVs.slice(i) = cholVsi;
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
                                    muJLsi += Bsjksi / Fs(j) * (AlphaJLsj - arma::kron(Bsj, EyeO) * AlphaJLNsj);
                                    //arma::kron(Bsj, EyeO) can be treated as arma::sp_mat as well but probably no need as O is very small and often equals 1
                                }
                                muJLsi = KappaInv * muJLsi;
                            }
                        }
                        if (i > 0) {
                            if (i < h) nnMaxNum = i - 1;
                            else nnMaxNum = h - 1;
                            arma::mat Bsi = Bs(arma::span(i, i), arma::span(0, nnMaxNum));
                            arma::uvec nni = nnInd(arma::span(i, i), arma::span(0, nnMaxNum)).t();
                            arma::mat AlphaJLNsi = arma::reshape(AlphaJLMat.cols(nni), (O * (nnMaxNum + 1)), 1);
                            muJLsi += arma::kron((1 / Fs(i)) * Bsi, KappaInv) * AlphaJLNsi;
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
                        if (truncated) { 
                            AlphaJLsiNew = arma::trans(rtmvnormRcpp(Vsi * muJLsi, Vsi, lowerVec, upperVec)); 
                        }
                        else { 
                            arma::mat Vsi = Vs.slice(i);
                            arma::mat cholVsi = cholVs.slice(i);
                            AlphaJLsiNew = rmvnormRcppNew(1, Vsi * muJLsi, cholVsi);
                        }
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
std::pair<paraVaryLj, listParaVaryLj> SampleXi(datobjVaryLj DatObj, paraVaryLj Para, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat YStarWide = DatObj.YStarWide;
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
    //arma::field<arma::mat> Weights = ListPara.Weights;
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
    return std::pair<paraVaryLj, listParaVaryLj>(Para, ListPara);
}



//Function to sample theta_jl's using a Gibbs sampler step---------------------------------------------------------------
std::pair<paraVaryLj, listParaVaryLj> SampleTheta(datobjVaryLj DatObj, paraVaryLj Para, listParaVaryLj ListPara) {
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int Nu = DatObj.Nu;
    arma::mat YStarWide = DatObj.YStarWide;
    arma::mat EyeNu = DatObj.EyeNu;
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
    return std::pair<paraVaryLj, listParaVaryLj>(Para, ListPara);
}



//Function to sample Latent U using a Gibbs sampler step-------------------------------------------------------------------
paraVaryLj SampleU(datobjVaryLj DatObj, paraVaryLj Para, listParaVaryLj ListPara) {
    //Set data objects
    int M = DatObj.M;
    int K = DatObj.K;
    int O = DatObj.O;

    //Set parameters
    arma::field<arma::mat> Weights = ListPara.Weights;
    arma::umat Xi = Para.Xi;
    // arma::colvec LjVec = Para.LjVec;

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
