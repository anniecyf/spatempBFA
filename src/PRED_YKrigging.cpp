#include <RcppArmadillo.h>
#include "PRED_predictions.h"
#include "MCMC_bfaSpatTemp.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube YKriggingTemp(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose) {
    datobjPREDtemp DatObj = ConvertDatObjPREDtemp(DatObj_List);
    paraPREDtempY Para = ConvertParaPREDtempY(Para_List);
    arma::Col<int> FamilyInd = DatObj.FamilyInd;
    int M = DatObj.M;
    int O = DatObj.O;
    int C = DatObj.C;
    int K = DatObj.K;
    int P = DatObj.P;
    int NNewVisits = DatObj.NNewVisits;
    arma::cube Trials = DatObj.Trials; //of dim M * O * NNewVisits if C = 0, and of dim M * C * NNewVisits otherwise
    arma::mat NewX = DatObj.NewX;
    arma::mat LambdaMat = Para.Lambda;
    arma::mat Sigma2Mat = Para.Sigma2;
    arma::mat BetaMat = Para.Beta;
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Y: 0%.. ";
    }
    arma::cube Out(M * O, NNewVisits, NKeep);
    arma::mat YMax(M, O, arma::fill::zeros), Lambda(M * O, K), OutTemp(M * O, NNewVisits);
    arma::mat BigPhi(K, NNewVisits), Mean(M * O, NNewVisits), MeanMat(M, O), PredMat(M, O, arma::fill::zeros);
    arma::colvec MeanVec(M), Pi(M), TrialsVec(M), Beta(P);
    arma::umat ProbitOnes;
    int FamilyType, count;
    for (arma::uword s = 0; s < NKeep; s++) {
        BigPhi = arma::reshape(EtaKrig.col(s), K, NNewVisits);
        Lambda = arma::reshape(LambdaMat.row(s), K, M * O).t();
        Beta = BetaMat.row(s).t();
        arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), M, (O - C));
        Mean = Lambda * BigPhi + arma::reshape(NewX * Beta, M * O, NNewVisits);
        for (arma::uword n = 0; n < NNewVisits; n++) {
            MeanMat = arma::reshape(Mean.col(n), M, O);
            arma::mat TrialsMat = Trials.slice(n); //of dim M*O if C=0, and of dim M*C otherwise
            count = 0;
            for (arma::uword f = 0; f < O; f++) {
                FamilyType = FamilyInd(f);
                if (FamilyType != 3) {
                    arma::vec SD = arma::sqrt(Sigma2.col(count));
                    MeanVec = MeanMat.col(f);
                    PredMat.col(f) = rnormVecRcpp(MeanVec, SD);
                    count++;
                    if (FamilyType == 1) { //Probit
                        PredMat = arma::max(PredMat, YMax);
                        ProbitOnes = find(PredMat > 0);
                        PredMat(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
                    }
                    else if (FamilyType == 2) { //Tobit
                        PredMat.col(f) = arma::max(PredMat.col(f), YMax.col(f));
                    }
                }
                //Categorical
                else if (FamilyType == 3) {
                    MeanVec = MeanMat.col(f);
                    TrialsVec = TrialsMat.col(f - count);
                    Pi = arma::exp(MeanVec) / (1 + arma::exp(MeanVec));
                    for (arma::uword i = 0; i < M; i++) PredMat(i, f) = rbinomRcpp(TrialsVec(i), Pi(i));
                }
            }
            OutTemp.col(n) = arma::vectorise(PredMat);
        }
        Out.slice(s) = OutTemp;
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
    //Return PPD
    return Out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube VAR1YKriggingTemp(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose) {
    VAR1datobjPREDtemp DatObj = VAR1ConvertDatObjPREDtemp(DatObj_List);
    paraPREDtempY Para = ConvertParaPREDtempY(Para_List);
    arma::Col<int> FamilyInd = DatObj.FamilyInd;
    int M = DatObj.M;
    int O = DatObj.O;
    int C = DatObj.C;
    int K = DatObj.K;
    int P = DatObj.P;
    int NNewTime = DatObj.NNewTime;
    arma::cube Trials = DatObj.Trials; //of dim M * O * NNewTime if C = 0, and of dim M * C * NNewTime otherwise
    arma::mat NewX = DatObj.NewX;
    arma::mat LambdaMat = Para.Lambda;
    arma::mat Sigma2Mat = Para.Sigma2;
    arma::mat BetaMat = Para.Beta;
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Y: 0%.. ";
    }
    arma::cube Out(M * O, NNewTime, NKeep);
    arma::mat YMax(M, O, arma::fill::zeros), Lambda(M * O, K), OutTemp(M * O, NNewTime);
    arma::mat BigPhi(K, NNewTime), Mean(M * O, NNewTime), MeanMat(M, O), PredMat(M, O, arma::fill::zeros);
    arma::colvec MeanVec(M), Pi(M), TrialsVec(M), Beta(P);
    arma::umat ProbitOnes;
    int FamilyType, count;
    for (arma::uword s = 0; s < NKeep; s++) {
        BigPhi = arma::reshape(EtaKrig.col(s), K, NNewTime);
        Lambda = arma::reshape(LambdaMat.row(s), K, M * O).t();
        Beta = BetaMat.row(s).t();
        arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), M, (O - C));
        Mean = Lambda * BigPhi + arma::reshape(NewX * Beta, M * O, NNewTime);
        for (arma::uword n = 0; n < NNewTime; n++) {
            MeanMat = arma::reshape(Mean.col(n), M, O);
            arma::mat TrialsMat = Trials.slice(n); //of dim M*O if C=0, and of dim M*C otherwise
            count = 0;
            for (arma::uword f = 0; f < O; f++) {
                FamilyType = FamilyInd(f);
                if (FamilyType != 3) {
                    arma::vec SD = arma::sqrt(Sigma2.col(count));
                    MeanVec = MeanMat.col(f);
                    PredMat.col(f) = rnormVecRcpp(MeanVec, SD);
                    count++;
                    if (FamilyType == 1) { //Probit
                        PredMat = arma::max(PredMat, YMax);
                        ProbitOnes = find(PredMat > 0);
                        PredMat(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
                    }
                    else if (FamilyType == 2) { //Tobit
                        PredMat.col(f) = arma::max(PredMat.col(f), YMax.col(f));
                    }
                }
                //Categorical
                else if (FamilyType == 3) {
                    MeanVec = MeanMat.col(f);
                    TrialsVec = TrialsMat.col(f - count);
                    Pi = arma::exp(MeanVec) / (1 + arma::exp(MeanVec));
                    for (arma::uword i = 0; i < M; i++) PredMat(i, f) = rbinomRcpp(TrialsVec(i), Pi(i));
                }
            }
            OutTemp.col(n) = arma::vectorise(PredMat);
        }
        Out.slice(s) = OutTemp;
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
    //Return PPD
    return Out;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube YKriggingSpat(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat LambdaKrig, int NKeep, bool Verbose) {
    //discuss here
    datobjPREDspat DatObj = ConvertDatObjPREDspatVaryLj(DatObj_List);
    paraPREDspatY Para = ConvertParaPREDspatY(Para_List);
    arma::Col<int> FamilyInd = DatObj.FamilyInd;
    int O = DatObj.O;
    int C = DatObj.C;
    int K = DatObj.K;
    int P = DatObj.P;
    int T = DatObj.Nu;
    int NNewLoc = DatObj.NNewLoc;
    arma::cube Trials = DatObj.Trials; //of dim T * O * NNewLoc if C = 0, and of dim T * C * NNewLoc otherwise
    arma::mat NewX = DatObj.NewX;
    double sigma2a = DatObj.sigma2hyparaA;
    double sigma2b = DatObj.sigma2hyparaB;
    arma::mat EtaMat = Para.Eta;
    arma::mat BetaMat = Para.Beta;
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Y: 0%.. ";
    }
    arma::mat BigPhi(K, T), Lambda(NNewLoc * O, K), OutSpat(T * O, NNewLoc);
    arma::cube Out(T * O, NNewLoc, NKeep);
    arma::mat YMax(T, O, arma::fill::zeros);
    arma::mat Mean(T * O, NNewLoc), MeanMat(T, O), PredMat(T, O, arma::fill::zeros);
    arma::colvec MeanVec(T), Pi(T), TrialsVec(T);
    arma::umat ProbitOnes;
    int FamilyType, count;
    arma::mat Sigma2(T, (O - C));
    for (arma::uword t = 0; t < T; t++) {
        for (arma::uword o = 0; o < (O - C); o++) {
            Sigma2(t, o) = rigammaRcpp(sigma2a, sigma2b);
        }
    }
    for (arma::uword s = 0; s < NKeep; s++) {
        BigPhi = arma::reshape(EtaMat.row(s), K, T);
        Lambda = arma::reshape(LambdaKrig.col(s), NNewLoc * O, K);//LambdaKrig is of dim (NNewLoc * O * K) x NKeep
        //rows of Lambda first ordered by observation type and then spatially, same as that for rows of NewX corresponding to any t
        //(the first O rows correspond to the first new location point for prediction, the next O rows correspond to the second new location point for prediction and so on)              
        if (P == 0) {
            Mean = arma::reshape(arma::trans(Lambda * BigPhi), T * O, NNewLoc);
            //rows of Mean first ordered temporally and then by observation type (the first T rows correspond to the first observation type, the next T rows correspond to the second observation type and so on)
        } 
        else {// if P > 0
            arma::colvec Beta = BetaMat.row(s).t();
            Mean = arma::reshape(arma::trans(Lambda * BigPhi + arma::reshape(NewX * Beta, NNewLoc * O, T)), T * O, NNewLoc); //Lambda * BigPhi + arma::reshape(NewX * Beta, NNewLoc * O, T) is of dim (NNewLoc * O) x T
        }        
        //Mean = Lambda * BigPhi + arma::reshape(NewX * Beta, NNewLoc * O, T);
        for (arma::uword n = 0; n < NNewLoc; n++) {
            MeanMat = arma::reshape(Mean.col(n), T, O);
            arma::mat TrialsMat = Trials.slice(n); //of dim T*O if C=0, and of dim T*C otherwise
            count = 0;
            for (arma::uword f = 0; f < O; f++) {
                FamilyType = FamilyInd(f);
                if (FamilyType != 3) {
                    arma::vec SD = arma::sqrt(Sigma2.col(count));//of length T
                    MeanVec = MeanMat.col(f);
                    PredMat.col(f) = rnormVecRcpp(MeanVec, SD);
                    count++;
                    if (FamilyType == 1) { //Probit
                        PredMat = arma::max(PredMat, YMax);
                        ProbitOnes = find(PredMat > 0);
                        PredMat(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
                    }
                    else if (FamilyType == 2) { //Tobit
                        PredMat.col(f) = arma::max(PredMat.col(f), YMax.col(f));
                    }
                }
                //Categorical
                else {//if FamilyType == 3
                    MeanVec = MeanMat.col(f);
                    TrialsVec = TrialsMat.col(f - count);
                    Pi = arma::exp(MeanVec) / (1 + arma::exp(MeanVec));
                    for (arma::uword t = 0; t < T; t++) PredMat(t, f) = rbinomRcpp(TrialsVec(t), Pi(t));
                }
            }
            OutSpat.col(n) = arma::vectorise(PredMat);//the first T predictions correspond to the first obs type, the next T predictions correspond to the second obs type and so on
        }
        Out.slice(s) = OutSpat;
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
    //Return PPD
    return Out;
}

