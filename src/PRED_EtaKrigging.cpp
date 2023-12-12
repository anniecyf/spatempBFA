#include <RcppArmadillo.h>
#include "PRED_predictions.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {
    datobjPREDtemp DatObj = ConvertDatObjPREDtemp(DatObj_List);
    paraPREDtempEta Para = ConvertParaPREDtempEta(Para_List);
    int Nu = DatObj.Nu;
    int K = DatObj.K;
    int ET = DatObj.ET;
    int IT = DatObj.IT;
    int NNewVisits = DatObj.NNewVisits;
    int TempCorInd = DatObj.TempCorInd;
    int seasonPeriod = DatObj.seasonPeriod;
    arma::mat TimeDist = DatObj.TimeDist;
    arma::uvec NewVisits = DatObj.NewVisits;
    arma::uvec OriginalVisits = DatObj.OriginalVisits;
    arma::mat EyeK(K, K, arma::fill::eye);
    arma::mat PsiMat = Para.Psi;
    arma::mat UpsilonMat = Para.Upsilon;
    arma::mat EtaMat = Para.Eta;
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Eta: 0%.. ";
    }
    arma::mat Upsilon(K, K), Sigma(K * NNewVisits, K * NNewVisits);
    arma::rowvec UpsilonRow(K);
    arma::colvec Eta(K * Nu), Mean(K * NNewVisits);
    arma::mat HPsiFull(Nu + NNewVisits, Nu + NNewVisits), KrigOut(K * NNewVisits, NKeep);
    arma::mat CovNewOrig(NNewVisits, Nu), CovOrigOrig(Nu, Nu), CovOrigOrigInv(Nu, Nu), CovPlus(NNewVisits, Nu);
    arma::mat CovNewNew(NNewVisits, NNewVisits), CovStar(NNewVisits, NNewVisits);
    double Psi;
    int counter;
    for (arma::uword s = 0; s < NKeep; s++) {
        Eta = EtaMat.row(s).t();
        UpsilonRow = UpsilonMat.row(s);
        counter = 0;
        for (arma::uword i = 0; i < K; i++) {
            for (arma::uword j = 0; j <= i; j++) {
                Upsilon(i, j) = UpsilonRow(counter);
                counter++;
            }
        }
        Upsilon = arma::symmatl(Upsilon);
        Psi = arma::as_scalar(PsiMat.row(s));
        if (IT == 0) {
            arma::mat EyeNew(NNewVisits, NNewVisits, arma::fill::eye);
            Sigma = arma::kron(EyeNew, Upsilon); //CovStar = EyeNew
            Mean = arma::colvec(NNewVisits * K, arma::fill::zeros); //CovPlus = matrix of zeros with proper dim
        }
        else if (ET == 1) { //i.e., when include.time = equalTimeDist = TRUE
            //Calculate covariance components
            HPsiFull = getH(Psi, TempCorInd, TimeDist, Nu + NNewVisits, seasonPeriod);
            CovNewOrig = HPsiFull(NewVisits, OriginalVisits);
            CovOrigOrigInv = getInvH(Psi, Nu, TempCorInd, seasonPeriod);
            CovPlus = CovNewOrig * CovOrigOrigInv;
            CovNewNew = HPsiFull(NewVisits, NewVisits);
            CovStar = CovNewNew - CovPlus * arma::trans(CovNewOrig);
            //Compute moments
            Sigma = arma::kron(CovStar, Upsilon);
            Mean = arma::kron(CovPlus, EyeK) * Eta;
        }
        else { //when include.time = TRUE and equalTimeDist = FALSE
            //Calculate covariance components
            HPsiFull = getH(Psi, TempCorInd, TimeDist, Nu + NNewVisits, seasonPeriod);
            CovNewOrig = HPsiFull(NewVisits, OriginalVisits);
            CovOrigOrig = HPsiFull(OriginalVisits, OriginalVisits);
            CovPlus = CovNewOrig * CholInv(CovOrigOrig);
            CovNewNew = HPsiFull(NewVisits, NewVisits);
            CovStar = CovNewNew - CovPlus * arma::trans(CovNewOrig);
            //Compute moments
            Sigma = arma::kron(CovStar, Upsilon);
            Mean = arma::kron(CovPlus, EyeK) * Eta;
        }
        KrigOut.col(s) = rmvnormRcpp(1, Mean, Sigma);
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
    return KrigOut;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat VAR1EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {
    VAR1datobjPREDtemp DatObj = VAR1ConvertDatObjPREDtemp(DatObj_List);
    VAR1paraPREDtempEta Para = VAR1ConvertParaPREDtempEta(Para_List);
    int Nu = DatObj.Nu;
    int K = DatObj.K;
    int NNewTime = DatObj.NNewTime;
    arma::mat Amat = Para.A; //NKeep x (K * K)
    arma::mat UpsilonMat = Para.Upsilon;
    arma::mat EtaMat = Para.Eta;
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Eta: 0%.. ";
    }
    arma::mat Upsilon(K, K), A(K, K);
    arma::rowvec UpsilonRow(K);
    arma::colvec Eta(K * Nu), etaPrev(K), etaPredColvec(K);
    arma::mat KrigOut(K * NNewTime, NKeep), krigEtaPerIter(K, NNewTime);
    int counter;
    for (arma::uword s = 0; s < NKeep; s++) {
        Eta = EtaMat.row(s).t();
        UpsilonRow = UpsilonMat.row(s);
        counter = 0;
        for (arma::uword i = 0; i < K; i++) {
            for (arma::uword j = 0; j <= i; j++) {
                Upsilon(i, j) = UpsilonRow(counter);
                counter++;
            }
        }
        Upsilon = arma::symmatl(Upsilon);
        A = arma::trans(arma::reshape(Amat.row(s), K, K));
        etaPrev = Eta(arma::span(K * (Nu - 1), (K * Nu - 1))); //initial value: etaT
        for (arma::uword t = 0; t < NNewTime; t++) {
            etaPredColvec = rmvnormRcpp(1, A * etaPrev, Upsilon);
            krigEtaPerIter.col(t) = etaPredColvec;
            etaPrev = etaPredColvec;
        }       
        KrigOut.col(s) = arma::vectorise(krigEtaPerIter);
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
    return KrigOut;
}
