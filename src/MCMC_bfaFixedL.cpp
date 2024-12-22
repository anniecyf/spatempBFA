#include <RcppArmadillo.h>
#include "MCMC_bfaSpatemp.h"
#include <ctime>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List bfaRcppFixedL(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List ParaCL_List, Rcpp::List SpatPara_List,
                       Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive) {
    //Convert Rcpp::Lists to C++ structs and set objects to be used in MCMC sampler 
    datobjFixedL DatObj = ConvertDatObjFixedL(DatObj_List);
    hypara HyPara = ConvertHyPara(HyPara_List);
    metrobj MetrObj = ConvertMetrObj(MetrObj_List);
    paraFixedL Para = ConvertParaFixedL(Para_List);
    int BNP = DatObj.CL;
    paraCLfixedL ParaCL {};
    if (BNP == 1) ParaCL = ConvertParaCLfixedL(ParaCL_List);
    dataug DatAug = ConvertDatAug(DatAug_List);
    mcmcobj McmcObj = ConvertMcmcObj(McmcObj_List);
    int IS = DatObj.IS;
    int SpCorInd = DatObj.SpCorInd;
    int scenario;
    spatpara1 SpatPara1 {};
    spatpara2 SpatPara2 {};
    spatpara3 SpatPara3 {};
    if (!(IS == 1 && SpCorInd == 0 && DatObj.SA == 1)) scenario = 1;
    else if (BNP == 0) scenario = 2;
    else if (DatObj.alphaMethodInd == 1) scenario = 3; //alphaMethod = "sequential"
    else scenario = 4; //alphaMethod = "block"
    if (scenario == 1) SpatPara1 = ConvertSpatPara1(SpatPara_List);
    else if (scenario == 2) SpatPara2 = ConvertSpatPara2(SpatPara_List);
    else if (scenario == 3) SpatPara3 = ConvertSpatPara3(SpatPara_List);
    else SpatPara2 = ConvertSpatPara2(SpatPara_List);
    arma::Col<int> FamilyInd = DatObj.FamilyInd;
    int NTotal = McmcObj.NTotal;
    int NBurn = McmcObj.NBurn;
    int NKeep = McmcObj.NKeep;
    arma::vec WhichPilotAdapt = McmcObj.WhichPilotAdapt;
    arma::vec WhichKeep = McmcObj.WhichKeep;
    arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
    arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
    arma::vec WhichSamplerProgress = McmcObj.WhichSamplerProgress;
    std::pair<paraFixedL, paraCLfixedL> XiThetaUpdate {};
    std::pair<datobjFixedL, paraFixedL> DAUpdate {};
    std::pair<paraFixedL, metrobj> Update {};
    std::pair<spatpara1, metrobj> Update1 {};
    std::pair<spatpara2, metrobj> Update2 {};
    std::pair<spatpara3, metrobj> Update3 {};
    arma::mat GibbsStepTime(10, NKeep, arma::fill::zeros);
    //std::time_t t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
    arma::colvec GibbsStepTimeVec(10, arma::fill::zeros);//auto time1 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    //ctime double time1 = std::difftime(t2, t1);// in seconds
    int keptIter = 0;
    bool storePostPara = false;

    //User output
    BeginBurnInProgress(McmcObj, Interactive);

    //Begin MCMC Sampler
    for (int s = 1; s < NTotal + 1; s++) {
        if (std::find(WhichKeep.begin(), WhichKeep.end(), s) != WhichKeep.end()) {
            keptIter += 1;
            storePostPara = true;
        }
        else {
            storePostPara = false;
        }
        if (s % 50 == 0) Rcpp::checkUserInterrupt();

        // Data Augmentation Step
        if (any(FamilyInd != 0)) {
            DAUpdate = SampleY(DatObj, Para, DatAug);
            DatObj = DAUpdate.first;
            Para = DAUpdate.second;
        }

        if (BNP == 1) {
            //Gibbs step for Xi
            auto t1 = std::chrono::high_resolution_clock::now();
            XiThetaUpdate = SampleXi(DatObj, Para, ParaCL);
            auto t2 = std::chrono::high_resolution_clock::now();
            Para = XiThetaUpdate.first;
            ParaCL = XiThetaUpdate.second;

            //Gibbs step for Theta
            auto t3 = std::chrono::high_resolution_clock::now();
            XiThetaUpdate = SampleTheta(DatObj, Para, ParaCL);
            auto t4 = std::chrono::high_resolution_clock::now();
            Para = XiThetaUpdate.first;
            ParaCL = XiThetaUpdate.second;

            //Gibbs step for Delta
            auto t5 = std::chrono::high_resolution_clock::now();
            ParaCL = SampleDelta(DatObj, Para, HyPara, ParaCL);

            //Gibbs step for Z
            auto t6 = std::chrono::high_resolution_clock::now();
            ParaCL = SampleZ(DatObj, Para, ParaCL);
            auto t7 = std::chrono::high_resolution_clock::now();

            GibbsStepTimeVec(0) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t7 - t6).count()); // for z
            GibbsStepTimeVec(1) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()); // for xi
            GibbsStepTimeVec(2) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count()); // for theta
            GibbsStepTimeVec(3) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count()); // for delta

            //Gibbs steps for Alpha and Kappa
            if (scenario == 1) {
                auto t8 = std::chrono::high_resolution_clock::now();
                ParaCL = SampleAlpha(DatObj, Para, ParaCL, SpatPara1, keptIter, storePostPara);
                auto t9 = std::chrono::high_resolution_clock::now();
                Para = SampleKappa(DatObj, Para, ParaCL, SpatPara1, HyPara);
                auto t10 = std::chrono::high_resolution_clock::now();
                GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
                GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
            }
            else if (scenario == 2 || scenario == 4) {
                auto t8 = std::chrono::high_resolution_clock::now();
                ParaCL = SampleAlpha(DatObj, Para, ParaCL, SpatPara2, keptIter, storePostPara);
                auto t9 = std::chrono::high_resolution_clock::now();
                Para = SampleKappa(DatObj, Para, ParaCL, SpatPara2, HyPara);
                auto t10 = std::chrono::high_resolution_clock::now();
                GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
                GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
            }
            else {//when scenario == 3 (when IS = 1, SpCorInd = 0, DatObj.SA = 1, BNP = 1, and alphaMethod = "sequential")
                auto t8 = std::chrono::high_resolution_clock::now();
                ParaCL = SampleAlpha(DatObj, Para, ParaCL, SpatPara3, keptIter, storePostPara);
                auto t9 = std::chrono::high_resolution_clock::now();
                Para = SampleKappa(DatObj, Para, ParaCL, SpatPara3, HyPara);
                auto t10 = std::chrono::high_resolution_clock::now();
                GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
                GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
            }
            //GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
            //GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
        }
        else {// when BNP == 0
            //Gibbs steps for Lambda and Kappa
            if (scenario == 1) {
                auto t8 = std::chrono::high_resolution_clock::now();
                Para = SampleLambda(DatObj, Para, SpatPara1);
                auto t9 = std::chrono::high_resolution_clock::now();
                Para = SampleKappa(DatObj, Para, SpatPara1, HyPara);
                auto t10 = std::chrono::high_resolution_clock::now();
                GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
                GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
            }
            else {//scenario must == 2 in this case
                auto t8 = std::chrono::high_resolution_clock::now();
                Para = SampleLambda(DatObj, Para, SpatPara2);
                auto t9 = std::chrono::high_resolution_clock::now();
                Para = SampleKappa(DatObj, Para, SpatPara2, HyPara);
                auto t10 = std::chrono::high_resolution_clock::now();
                GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for alpha
                GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
            }
            //GibbsStepTimeVec(4) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t9 - t8).count()); // for lambda
            //GibbsStepTimeVec(5) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t10 - t9).count()); // for kappa
        }

        //Metropolis step for Rho
        if (IS == 1) {
            if (scenario == 1) {// when spatApprox = FALSE; need to discuss whether spatial.structure is "discrete" or "continuous"
                auto t11 = std::chrono::high_resolution_clock::now();
                Update1 = SampleRho(DatObj, Para, ParaCL, SpatPara1, HyPara, MetrObj);
                auto t12 = std::chrono::high_resolution_clock::now();
                SpatPara1 = Update1.first;
                MetrObj = Update1.second;
                GibbsStepTimeVec(6) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count()); // for rho
            }
            else if (scenario == 2 || scenario == 4) {
                auto t11 = std::chrono::high_resolution_clock::now();
                Update2 = SampleRho(DatObj, Para, ParaCL, SpatPara2, HyPara, MetrObj);
                auto t12 = std::chrono::high_resolution_clock::now();
                SpatPara2 = Update2.first;
                MetrObj = Update2.second;
                GibbsStepTimeVec(6) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t12 - t11).count()); // for rho
            }
            else {//when scenario == 3 (when IS = 1, SpCorInd = 0, DatObj.SA = 1, BNP = 1, and alphaMethod = "sequential")
                auto t11 = std::chrono::high_resolution_clock::now();
                Update3 = SampleRho(DatObj, Para, ParaCL, SpatPara3, HyPara, MetrObj);
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

        //Gibbs step for Upsilon
        auto t14 = std::chrono::high_resolution_clock::now();
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
        GibbsStepTimeVec(7) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t14 - t13).count()); // for eta
        GibbsStepTimeVec(8) = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t15 - t14).count()); // for upsilon       
        int numTotalPara = size(RawSamples, 0);
        if (storePostPara == true) {
            int keepInd = arma::as_scalar(arma::find(s == WhichKeep, 1));
            GibbsStepTime.col(keepInd) = GibbsStepTimeVec;
            if (BNP == 0) {
                if (scenario == 1) RawSamples.col(keepInd) = StoreSamples(DatObj, Para, SpatPara1, numTotalPara);
                else if (scenario == 2 || scenario == 4) RawSamples.col(keepInd) = StoreSamples(DatObj, Para, SpatPara2, numTotalPara);
                else RawSamples.col(keepInd) = StoreSamples(DatObj, Para, SpatPara3, numTotalPara);
            }
            else {
                if (scenario == 1) RawSamples.col(keepInd) = StoreSamples(DatObj, Para, ParaCL, SpatPara1, numTotalPara);
                else if (scenario == 2 || scenario == 4) RawSamples.col(keepInd) = StoreSamples(DatObj, Para, ParaCL, SpatPara2, numTotalPara);
                else RawSamples.col(keepInd) = StoreSamples(DatObj, Para, ParaCL, SpatPara3, numTotalPara);
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
    
    //Return raw samples
    return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                              Rcpp::Named("GibbsStepTime") = GibbsStepTime,
                              Rcpp::Named("metropolis") = Metropolis);

    //End MCMC sampler function  
}
