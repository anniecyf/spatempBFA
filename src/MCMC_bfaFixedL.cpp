#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"

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
            XiThetaUpdate = SampleXi(DatObj, Para, ParaCL);
            Para = XiThetaUpdate.first;
            ParaCL = XiThetaUpdate.second;

            //Gibbs step for Theta
            XiThetaUpdate = SampleTheta(DatObj, Para, ParaCL);
            Para = XiThetaUpdate.first;
            ParaCL = XiThetaUpdate.second;

            //Gibbs step for Delta
            ParaCL = SampleDelta(DatObj, Para, HyPara, ParaCL);

            //Gibbs step for Z
            ParaCL = SampleZ(DatObj, Para, ParaCL);

            //Gibbs steps for Alpha and Kappa
            if (scenario == 1) {
                ParaCL = SampleAlpha(DatObj, Para, ParaCL, SpatPara1, keptIter, storePostPara);
                Para = SampleKappa(DatObj, Para, ParaCL, SpatPara1, HyPara);
            }
            else if (scenario == 2 || scenario == 4) {
                ParaCL = SampleAlpha(DatObj, Para, ParaCL, SpatPara2, keptIter, storePostPara);
                Para = SampleKappa(DatObj, Para, ParaCL, SpatPara2, HyPara);
            }
            else {//when scenario == 3 (when IS = 1, SpCorInd = 0, DatObj.SA = 1, BNP = 1, and alphaMethod = "sequential")
                ParaCL = SampleAlpha(DatObj, Para, ParaCL, SpatPara3, keptIter, storePostPara);
                Para = SampleKappa(DatObj, Para, ParaCL, SpatPara3, HyPara);
            }
        }
        else {// when BNP == 0
            //Gibbs steps for Lambda and Kappa
            if (scenario == 1) {
                Para = SampleLambda(DatObj, Para, SpatPara1);
                Para = SampleKappa(DatObj, Para, SpatPara1, HyPara);
            }
            else {//scenario must == 2 in this case
                Para = SampleLambda(DatObj, Para, SpatPara2);
                Para = SampleKappa(DatObj, Para, SpatPara2, HyPara);
            }
        }

        //Metropolis step for Rho
        if (IS == 1) {
            if (scenario == 1) {// when spatApprox = FALSE; need to discuss whether spatial.structure is "discrete" or "continuous"
                Update1 = SampleRho(DatObj, Para, ParaCL, SpatPara1, HyPara, MetrObj);
                SpatPara1 = Update1.first;
                MetrObj = Update1.second;
            }
            else if (scenario == 2 || scenario == 4) {
                Update2 = SampleRho(DatObj, Para, ParaCL, SpatPara2, HyPara, MetrObj);
                SpatPara2 = Update2.first;
                MetrObj = Update2.second;
            }
            else {//when scenario == 3 (when IS = 1, SpCorInd = 0, DatObj.SA = 1, BNP = 1, and alphaMethod = "sequential")
                Update3 = SampleRho(DatObj, Para, ParaCL, SpatPara3, HyPara, MetrObj);
                SpatPara3 = Update3.first;
                MetrObj = Update3.second;
            }
        }

        //Gibbs step for Eta
        Para = SampleEta(DatObj, Para, HyPara);

        //Gibbs step for Upsilon
        Para = SampleUpsilon(DatObj, Para, HyPara);

        // Metropolis step for Psi
        if (DatObj.IT == 1) {
            Update = SamplePsi(DatObj, Para, HyPara, MetrObj);
            Para = Update.first;
            MetrObj = Update.second;
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

        //Store raw samples
        int numTotalPara = size(RawSamples, 0);
        if (storePostPara == true) {
            if (BNP == 0) {
                if (scenario == 1) RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para, SpatPara1, numTotalPara);
                else if (scenario == 2 || scenario == 4) RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para, SpatPara2, numTotalPara);
                else RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para, SpatPara3, numTotalPara);
            }
            else {
                if (scenario == 1) RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para, ParaCL, SpatPara1, numTotalPara);
                else if (scenario == 2 || scenario == 4) RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para, ParaCL, SpatPara2, numTotalPara);
                else RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para, ParaCL, SpatPara3, numTotalPara);
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
                              Rcpp::Named("metropolis") = Metropolis);

    //End MCMC sampler function  
}
