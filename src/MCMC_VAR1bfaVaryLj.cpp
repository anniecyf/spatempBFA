#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List VAR1bfaRcppVaryingLjs(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List SpatPara_List,
                       Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive) {

  //Convert Rcpp::Lists to C++ structs and set objects to be used in MCMC sampler 
  VAR1datobjVaryLj DatObj = VAR1ConvertDatObjVaryLj(DatObj_List);
  VAR1hypara HyPara = VAR1ConvertHyPara(HyPara_List);
  VAR1metrobj MetrObj = VAR1ConvertMetrObj(MetrObj_List);
  VAR1paraVaryLj Para = VAR1ConvertParaVaryLj(Para_List);
  dataug DatAug = ConvertDatAug(DatAug_List);
  mcmcobj McmcObj = ConvertMcmcObj(McmcObj_List);
  std::pair<listParaVaryLj, VAR1paraVaryLj> initUpdate = InitializeListPara(DatObj, Para);
  listParaVaryLj ListPara = initUpdate.first;
  Para = initUpdate.second;
  int IS = DatObj.IS;
  int SpCorInd = DatObj.SpCorInd;
  int K = DatObj.K;
  int alphaSequenInd = DatObj.alphaSequenInd;
  int scenario;
  spatpara1VaryLj SpatPara1;
  spatpara2VaryLj SpatPara2;
  spatpara3 SpatPara3;
  if (!(IS == 1 && SpCorInd == 0 && DatObj.SA == 1)) scenario = 1;
  else if (alphaSequenInd == 0) scenario = 2;
  else scenario = 3;
  if (scenario == 1) SpatPara1 = ConvertSpatPara1VaryLj(SpatPara_List);
  else if (scenario == 2) SpatPara2 = ConvertSpatPara2VaryLj(SpatPara_List);
  else SpatPara3 = ConvertSpatPara3(SpatPara_List);
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int spatPred = DatObj.spatPred;
  int storeW = DatObj.storeW;
  int NTotal = McmcObj.NTotal;
  int NBurn = McmcObj.NBurn;
  int NKeep = McmcObj.NKeep;
  arma::vec WhichPilotAdapt = McmcObj.WhichPilotAdapt;
  arma::vec WhichKeep = McmcObj.WhichKeep;
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress = McmcObj.WhichSamplerProgress;
  std::pair<VAR1datobjVaryLj, VAR1paraVaryLj> DAUpdate{};
  std::pair<spatpara1VaryLj, VAR1metrobj> Update1{};
  std::pair<spatpara2VaryLj, VAR1metrobj> Update2{};
  std::pair<spatpara3, VAR1metrobj> Update3{};
  std::pair<VAR1paraVaryLj, listParaVaryLj> ParaTogUpdate{}; 
  arma::field<arma::vec> ThetaSamples(K, NKeep);
  arma::field<arma::mat> AlphaSamples(K, NKeep);
  arma::field<arma::mat> WeightsSamples(K, NKeep);
  
  //User output
  BeginBurnInProgress(McmcObj, Interactive);

  //Begin MCMC Sampler
  for (int s = 1; s < NTotal + 1; s++) {

    //Check for user interrupt every 50 iterations
    if (s % 50 == 0) Rcpp::checkUserInterrupt();
    
    // Data Augmentation Step
    if (any(FamilyInd != 0)) {
      DAUpdate = SampleY(DatObj, Para, DatAug);
      DatObj = DAUpdate.first;
      Para = DAUpdate.second;
    }
    
    Para = SampleU(DatObj, Para, ListPara);

    //Gibbs step for Xi and update LjVec
    ParaTogUpdate = SampleXi(DatObj, Para, ListPara);
    Para = ParaTogUpdate.first;
    ListPara = ParaTogUpdate.second;

    //Gibbs step for Theta
    ParaTogUpdate = SampleTheta(DatObj, Para, ListPara);
    Para = ParaTogUpdate.first;
    ListPara = ParaTogUpdate.second;

    //Gibbs step for Delta
    Para = SampleDelta(DatObj, Para, HyPara, ListPara);
        
    //Gibbs steps for Alpha (and calculate the new Weights) and Kappa
    if (scenario == 1) {
        ListPara = SampleAlpha(DatObj, Para, SpatPara1, ListPara);
        Para = SampleKappa(DatObj, Para, SpatPara1, HyPara, ListPara);
    }
    else if (scenario == 2) {
        ListPara = SampleAlpha(DatObj, Para, SpatPara2, ListPara);
        Para = SampleKappa(DatObj, Para, SpatPara2, HyPara, ListPara);
    }
    else {
        ListPara = SampleAlpha(DatObj, Para, SpatPara3, ListPara);
        Para = SampleKappa(DatObj, Para, SpatPara3, HyPara, ListPara);
    }
   
    //Metropolis step for Rho
    if (IS == 1) {   
        if (scenario == 1) {// when spatApprox = FALSE; need to discuss whether spatial.structure is "discrete" or "continuous"
            Update1 = SampleRho(DatObj, Para, SpatPara1, HyPara, MetrObj, ListPara);
            SpatPara1 = Update1.first;
            MetrObj = Update1.second;
        }
        else if (scenario == 2) {
            Update2 = SampleRho(DatObj, Para, SpatPara2, HyPara, MetrObj, ListPara);
            SpatPara2 = Update2.first;
            MetrObj = Update2.second;
        }
        else {
            Update3 = SampleRho(DatObj, Para, SpatPara3, HyPara, MetrObj, ListPara);
            SpatPara3 = Update3.first;
            MetrObj = Update3.second;
        }
    }
    
    //Gibbs step for Eta
    Para = SampleEta(DatObj, Para, HyPara);
    
    //Gibbs step for Upsilon
    Para = SampleUpsilon(DatObj, Para, HyPara);
    
    // Metropolis step for A
    Para = SampleA(DatObj, Para, HyPara);

    //Gibbs step for Beta
    if (DatObj.P > 0) {
        Para = SampleBeta(DatObj, Para, HyPara);
    }
     
    //Gibbs step for Sigma2
    if (any(FamilyInd != 3)) {
      Para = SampleSigma2(DatObj, Para, HyPara);
    }
        
    //Pilot adaptation
    if (std::find(WhichPilotAdapt.begin(), WhichPilotAdapt.end(), s) != WhichPilotAdapt.end())
      MetrObj = PilotAdaptation(MetrObj, DatObj, McmcObj, s);
   
    //Store raw samples
    int numTotalPara = size(RawSamples, 0);
    if (std::find(WhichKeep.begin(), WhichKeep.end(), s) != WhichKeep.end()) {
        int keepInd = arma::as_scalar(arma::find(s == WhichKeep, 1));
        if (scenario == 1) { RawSamples.col(keepInd) = StoreSamples(DatObj, Para, SpatPara1, numTotalPara); }
        else if (scenario == 2) { RawSamples.col(keepInd) = StoreSamples(DatObj, Para, SpatPara2, numTotalPara); }
        else { RawSamples.col(keepInd) = StoreSamples(DatObj, Para, SpatPara3, numTotalPara); }
        if (spatPred == 1) {
            ThetaSamples.col(keepInd) = ListPara.Theta;
            AlphaSamples.col(keepInd) = ListPara.Alpha;
        }
        if (storeW == 1) {
            WeightsSamples.col(keepInd) = ListPara.Weights;            
        }     
    }
        
    //Update burn-in progress bar
    if (Interactive) if (std::find(WhichBurnInProgress.begin(), WhichBurnInProgress.end(), s) != WhichBurnInProgress.end())
      UpdateBurnInBar(s, McmcObj);
    if (!Interactive) if (std::find(WhichBurnInProgressInt.begin(), WhichBurnInProgressInt.end(), s) != WhichBurnInProgressInt.end())
      UpdateBurnInBarInt(s, McmcObj);

    //Post burn-in progress
    if (s == NBurn) Rcpp::Rcout << std::fixed << "\nSampler progress:  0%.. ";
    if (std::find(WhichSamplerProgress.begin(), WhichSamplerProgress.end(), s) != WhichSamplerProgress.end())
       SamplerProgress(s, McmcObj);

  //End MCMC Sampler
  }

  //Output Metropolis object for summary
  Rcpp::List Metropolis = OutputMetrObj(MetrObj);

  //Return raw samples alphafield
  if (spatPred == 1) {
      if (storeW == 1) {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("metropolis") = Metropolis,
                                    Rcpp::Named("thetafield") = ThetaSamples,
                                    Rcpp::Named("alphafield") = AlphaSamples,
                                    Rcpp::Named("weightsfield") = WeightsSamples);
      }
      else {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("metropolis") = Metropolis,
                                    Rcpp::Named("thetafield") = ThetaSamples,
                                    Rcpp::Named("alphafield") = AlphaSamples);
      }
  }
  else { // when spatPara = 0 
      if (storeW == 1) {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("metropolis") = Metropolis,
                                    Rcpp::Named("weightsfield") = WeightsSamples);
      }
      else {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("metropolis") = Metropolis);
      }
  }

//End MCMC sampler function
}
