#include <RcppArmadillo.h>
#include "MCMC_bfaSpatemp.h"
#include <ctime>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List bfaRcppVaryingLjs(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List SpatPara_List,
                       Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive) {

  //Convert Rcpp::Lists to C++ structs and set objects to be used in MCMC sampler 
  datobjVaryLj DatObj = ConvertDatObjVaryLj(DatObj_List);
  hypara HyPara = ConvertHyPara(HyPara_List);
  metrobj MetrObj = ConvertMetrObj(MetrObj_List);
  paraVaryLj Para = ConvertParaVaryLj(Para_List);
  dataug DatAug = ConvertDatAug(DatAug_List);
  mcmcobj McmcObj = ConvertMcmcObj(McmcObj_List);
  std::pair<listParaVaryLj, paraVaryLj> initUpdate = InitializeListPara(DatObj, Para);
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
  std::pair<datobjVaryLj, paraVaryLj> DAUpdate{};
  std::pair<paraVaryLj, metrobj> Update{};
  std::pair<spatpara1VaryLj, metrobj> Update1{};
  std::pair<spatpara2VaryLj, metrobj> Update2{};
  std::pair<spatpara3, metrobj> Update3{};
  std::pair<paraVaryLj, listParaVaryLj> ParaTogUpdate{}; 
  arma::field<arma::vec> ThetaSamples(K, NKeep);
  arma::field<arma::mat> AlphaSamples(K, NKeep);
  arma::field<arma::mat> WeightsSamples(K, NKeep);
  arma::mat GibbsStepTime(10, NKeep, arma::fill::zeros);
  //std::time_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  arma::colvec GibbsStepTimeVec(10, arma::fill::zeros);//auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  //ctime double time1 = std::difftime(t2, t1);// in seconds
   
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
    
    auto t1 = std::chrono::high_resolution_clock::now();
    Para = SampleU(DatObj, Para, ListPara);
    auto t2 = std::chrono::high_resolution_clock::now();

    //Gibbs step for Xi and update LjVec
    ParaTogUpdate = SampleXi(DatObj, Para, ListPara);
    auto t3 = std::chrono::high_resolution_clock::now();
    Para = ParaTogUpdate.first;
    ListPara = ParaTogUpdate.second;

    //Gibbs step for Theta
    auto t4 = std::chrono::high_resolution_clock::now();
    ParaTogUpdate = SampleTheta(DatObj, Para, ListPara);
    auto t5 = std::chrono::high_resolution_clock::now();
    Para = ParaTogUpdate.first;
    ListPara = ParaTogUpdate.second;

    //Gibbs step for Delta
    auto t6 = std::chrono::high_resolution_clock::now();
    Para = SampleDelta(DatObj, Para, HyPara, ListPara);
    auto t7 = std::chrono::high_resolution_clock::now();
        
    //Gibbs steps for Alpha (and calculate the new Weights) and Kappa
    if (scenario == 1) {
        auto t8 = std::chrono::high_resolution_clock::now();
        ListPara = SampleAlpha(DatObj, Para, SpatPara1, ListPara);
        auto t9 = std::chrono::high_resolution_clock::now();
        Para = SampleKappa(DatObj, Para, SpatPara1, HyPara, ListPara);
        auto t10 = std::chrono::high_resolution_clock::now();
        GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
        GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
    }
    else if (scenario == 2) {
        auto t8 = std::chrono::high_resolution_clock::now();
        ListPara = SampleAlpha(DatObj, Para, SpatPara2, ListPara);
        auto t9 = std::chrono::high_resolution_clock::now();
        Para = SampleKappa(DatObj, Para, SpatPara2, HyPara, ListPara);
        auto t10 = std::chrono::high_resolution_clock::now();
        GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
        GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
    }
    else {
        auto t8 = std::chrono::high_resolution_clock::now();
        ListPara = SampleAlpha(DatObj, Para, SpatPara3, ListPara);
        auto t9 = std::chrono::high_resolution_clock::now();
        Para = SampleKappa(DatObj, Para, SpatPara3, HyPara, ListPara);
        auto t10 = std::chrono::high_resolution_clock::now();
        GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
        GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
    }
   
    //Metropolis step for Rho
    if (IS == 1) {   
        if (scenario == 1) {// when spatApprox = FALSE; need to discuss whether spatial.structure is "discrete" or "continuous"
            auto t11 = std::chrono::high_resolution_clock::now();
            Update1 = SampleRho(DatObj, Para, SpatPara1, HyPara, MetrObj, ListPara);
            auto t12 = std::chrono::high_resolution_clock::now();
            SpatPara1 = Update1.first;
            MetrObj = Update1.second;
            GibbsStepTimeVec(6) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count()); // for rho
        }
        else if (scenario == 2) {
            auto t11 = std::chrono::high_resolution_clock::now();
            Update2 = SampleRho(DatObj, Para, SpatPara2, HyPara, MetrObj, ListPara);
            auto t12 = std::chrono::high_resolution_clock::now();
            SpatPara2 = Update2.first;
            MetrObj = Update2.second;
            GibbsStepTimeVec(6) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count()); // for rho
        }
        else {
            auto t11 = std::chrono::high_resolution_clock::now();
            Update3 = SampleRho(DatObj, Para, SpatPara3, HyPara, MetrObj, ListPara);
            auto t12 = std::chrono::high_resolution_clock::now();
            SpatPara3 = Update3.first;
            MetrObj = Update3.second;
            GibbsStepTimeVec(6) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count()); // for rho
        }
        //GibbsStepTimeVec(6) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count()); // for rho
    }
    
    //Gibbs step for Eta
    auto t13 = std::chrono::high_resolution_clock::now();
    Para = SampleEta(DatObj, Para, HyPara);
    auto t14 = std::chrono::high_resolution_clock::now();
    
    //Gibbs step for Upsilon
    Para = SampleUpsilon(DatObj, Para, HyPara);
    auto t15 = std::chrono::high_resolution_clock::now();
    
    // Metropolis step for Psi
    if (DatObj.IT == 1) {
        auto t16 = std::chrono::high_resolution_clock::now();
        Update = SamplePsi(DatObj, Para, HyPara, MetrObj);
        auto t17 = std::chrono::high_resolution_clock::now();
        Para = Update.first;
        MetrObj = Update.second;
        GibbsStepTimeVec(9) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t17 - t16).count()); // for psi
    }

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
   
    //Store raw samples and posterior sampling time for each Gibbs sampler step (except for beta and sigma2) at each kept post-burn-in MCMC iteration
    GibbsStepTimeVec(0) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()); // for u
    GibbsStepTimeVec(1) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()); // for xi
    GibbsStepTimeVec(2) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t5 - t4).count()); // for theta
    GibbsStepTimeVec(3) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t7 - t6).count()); // for delta
    //GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
    //GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
    GibbsStepTimeVec(7) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t14 - t13).count()); // for eta
    GibbsStepTimeVec(8) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t15 - t14).count()); // for upsilon
    int numTotalPara = size(RawSamples, 0);
    if (std::find(WhichKeep.begin(), WhichKeep.end(), s) != WhichKeep.end()) {     
        int keepInd = arma::as_scalar(arma::find(s == WhichKeep, 1));
        GibbsStepTime.col(keepInd) = GibbsStepTimeVec;
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
                                    Rcpp::Named("GibbsStepTime") = GibbsStepTime,
                                    Rcpp::Named("metropolis") = Metropolis,
                                    Rcpp::Named("thetafield") = ThetaSamples,
                                    Rcpp::Named("alphafield") = AlphaSamples,
                                    Rcpp::Named("weightsfield") = WeightsSamples);
      }
      else {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("GibbsStepTime") = GibbsStepTime,
                                    Rcpp::Named("metropolis") = Metropolis,
                                    Rcpp::Named("thetafield") = ThetaSamples,
                                    Rcpp::Named("alphafield") = AlphaSamples);
      }
  }
  else { // when spatPara = 0 
      if (storeW == 1) {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("GibbsStepTime") = GibbsStepTime,
                                    Rcpp::Named("metropolis") = Metropolis,
                                    Rcpp::Named("weightsfield") = WeightsSamples);
      }
      else {
          return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                                    Rcpp::Named("GibbsStepTime") = GibbsStepTime,
                                    Rcpp::Named("metropolis") = Metropolis);
      }
  }

//End MCMC sampler function
}
