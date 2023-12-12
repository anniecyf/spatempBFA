#include <iostream>
#include <fstream>
#include <string>
#include <RcppArmadillo.h>
#include "PRED_predictions.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat LambdaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {
    datobjPREDspat DatObj = ConvertDatObjPREDspatFixedL(DatObj_List);
    paraPREDspatFixedLnoCLlambda Para = ConvertParaPREDspatLambdaFixedLnoCL(Para_List);
    int T = DatObj.Nu;
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int IS = DatObj.IS;
    int SA = DatObj.SA;
    int h = DatObj.h;
    int NNewLoc = DatObj.NNewLoc;
    arma::umat nnIndpred = DatObj.nnIndpred; // of dim NNewLoc x h
    arma::mat dist = DatObj.SpDist;
    arma::mat distOrigNew = DatObj.distOrigNew;
    arma::mat distNewNew = DatObj.distNewNew;
    // if IS = 1, then SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
    arma::mat RhoMat = Para.Rho;
    arma::mat KappaMat = Para.Kappa;
    arma::mat LambdaMat = Para.Lambda; //of dim NKeep x (M x O x K)    
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Lambda: 0%.. ";
    }
    double rho;   
    arma::mat LambdaKrig(NNewLoc * O * K, NKeep); 
    //For each kept MCMC iteration, the corresponding predicted entries for the factor loadings matrix are ordered first by observation type, then spatially, and finally by factor. 
    //(the first (NNewLoc x O) entries correspond to factor 1, the next (NNewLoc x O) entries correspond to factor 2 and so on;
    //the first O rows correspond to the first new location point for prediction, the next O rows correspond to the second new location point for prediction and so on)  
    arma::mat LambdaOrig(M * O, K);
    arma::mat LambdaKrigMat(NNewLoc * O, K);//rows of LambdaKrigMat first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on), different as that for LambdaOrig 
    if (O == 1) {
        double kappa;
        for (arma::uword s = 0; s < NKeep; s++) {
            rho = arma::as_scalar(RhoMat.row(s));
            kappa = arma::as_scalar(KappaMat.row(s));
            LambdaOrig = arma::trans(arma::reshape(LambdaMat.row(s), K, M)); 
            if (IS == 0) {
                arma::colvec Mean(NNewLoc, arma::fill::zeros);
                arma::mat FrhoNewNew(NNewLoc, NNewLoc, arma::fill::eye);
                arma::mat Sigma = kappa * FrhoNewNew;
                LambdaKrigMat = rmvnormRcpp(K, Mean, Sigma);
            }
            else if (SA == 0) {// in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                //Calculate covariance components
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                arma::mat FrhoPlus = arma::trans(FrhoOrigNew) * CholInv(FrhoOrigOrig);
                arma::mat Sigma = kappa * (FrhoNewNew - FrhoPlus * FrhoOrigNew);
                for (arma::uword j = 0; j < K; j++) {
                    arma::colvec LambdaOrigJ = LambdaOrig.col(j);
                    arma::colvec MeanJ = FrhoPlus * LambdaOrigJ;
                    LambdaKrigMat.col(j) = rmvnormRcpp(1, MeanJ, Sigma);
                }
            }
            else { //when include.space = spatApprox = TRUE; in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::vec sigma2vec(NNewLoc);
                arma::mat Bmat(NNewLoc, h);
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    arma::uvec nni = arma::trans(nnIndpred.row(i));
                    arma::colvec FrhoOrigsi = FrhoOrigNew.col(i);
                    arma::colvec Frhosinsi = FrhoOrigsi(nni);//colvec
                    arma::mat Bsi = arma::trans(Frhosinsi) * CholInv(FrhoOrigOrig(nni, nni));
                    Bmat.row(i) = Bsi;
                    double Fsi = arma::as_scalar(1 - Bsi * Frhosinsi);
                    if (Fsi <= 0) Fsi = 0.00001;
                    sigma2vec(i) = Fsi * kappa;                   
                }
                for (arma::uword j = 0; j < K; j++) {
                    arma::colvec LambdaOrigJ = LambdaOrig.col(j);//of length M
                    arma::colvec LambdaKrigJ(NNewLoc);
                    for (arma::uword i = 0; i < NNewLoc; i++) {
                        arma::uvec nni = arma::trans(nnIndpred.row(i));
                        arma::colvec LambdaOrigJNsi = LambdaOrigJ(nni);
                        arma::mat Bsi = Bmat.row(i); // of dim 1 x h
                        double sdsi = sqrt(sigma2vec(i)); 
                        double meanJsi = arma::as_scalar(Bsi * LambdaOrigJNsi);
                        LambdaKrigJ(i) = arma::as_scalar(rnormRcpp(1, meanJsi, sdsi));
                    }
                    LambdaKrigMat.col(j) = LambdaKrigJ;
                }
            }
            LambdaKrig.col(s) = arma::reshape(LambdaKrigMat, NNewLoc * K, 1);
            if (Verbose) {
                Rcpp::Rcout.precision(0);
                if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                    Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
            }
        }
    }
    else { // if O > 1
        arma::mat EyeO(O, O, arma::fill::eye);
        arma::rowvec KappaRow(O * (O + 1) / 2);
        arma::mat Kappa(O, O);
        int counter;
        for (arma::uword s = 0; s < NKeep; s++) {
            rho = arma::as_scalar(RhoMat.row(s));
            KappaRow = KappaMat.row(s);
            counter = 0;
            for (arma::uword i = 0; i < O; i++) {
                for (arma::uword j = 0; j <= i; j++) {
                    Kappa(i, j) = KappaRow(counter);
                    counter++;
                }
            }
            Kappa = arma::symmatl(Kappa);
            LambdaOrig = arma::trans(arma::reshape(LambdaMat.row(s), K, M * O)); //rows of LambdaOrig here ordered first spatially and then by observation type 
            //(the first M corresponds to obs type 1, the next M corresponds to obs type 2 and so on); hence we need to reshape each column of LambdaOrig later
            if (IS == 0) {
                arma::colvec Mean(NNewLoc * O, arma::fill::zeros);
                arma::mat FrhoNewNew(NNewLoc, NNewLoc, arma::fill::eye);
                arma::mat Sigma = arma::kron(FrhoNewNew, Kappa);
                LambdaKrigMat = rmvnormRcpp(K, Mean, Sigma);
            }
            else if (SA == 0) {// in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                //Calculate covariance components
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                arma::mat FrhoPlus = arma::trans(FrhoOrigNew) * CholInv(FrhoOrigOrig);
                arma::mat Sigma = arma::kron((FrhoNewNew - FrhoPlus * FrhoOrigNew), Kappa);
                arma::mat FrhoPlusKron = arma::kron(FrhoPlus, EyeO);
                for (arma::uword j = 0; j < K; j++) {
                    arma::colvec LambdaOrigJ = LambdaOrig.col(j);
                    arma::mat LambdaOrigJMat = arma::trans(arma::reshape(LambdaOrigJ, M, O));
                    LambdaOrigJ = arma::reshape(LambdaOrigJMat, M * O, 1);
                    arma::colvec MeanJ = FrhoPlusKron * LambdaOrigJ;
                    LambdaKrigMat.col(j) = rmvnormRcpp(1, MeanJ, Sigma);
                }
            }
            else { //when include.space = spatApprox = TRUE; in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::cube CovCube(O, O, NNewLoc);
                arma::cube BCube(O, h * O, NNewLoc);
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    arma::uvec nni = arma::trans(nnIndpred.row(i));
                    arma::colvec FrhoOrigsi = FrhoOrigNew.col(i);
                    arma::mat Frhosinsi = FrhoOrigsi(nni);//colvec
                    arma::mat Bsi = arma::trans(Frhosinsi) * CholInv(FrhoOrigOrig(nni, nni));
                    double Fsi = arma::as_scalar(1 - Bsi * Frhosinsi);
                    if (Fsi <= 0) Fsi = 0.00001;
                    CovCube.slice(i) = Fsi * Kappa;
                    BCube.slice(i) = arma::kron(Bsi, EyeO);
                }
                for (arma::uword j = 0; j < K; j++) {
                    arma::colvec LambdaOrigJ = LambdaOrig.col(j);
                    arma::mat LambdaOrigJMat = arma::trans(arma::reshape(LambdaOrigJ, M, O)); // of dim O x M
                    arma::mat LambdaKrigJMat(O, NNewLoc);
                    for (arma::uword i = 0; i < NNewLoc; i++) {
                        arma::uvec nni = arma::trans(nnIndpred.row(i));
                        arma::mat LambdaOrigJNsiMat = LambdaOrigJMat.cols(nni); // of dim O x h
                        arma::colvec LambdaOrigJNsi = arma::reshape(LambdaOrigJNsiMat, O * h, 1);
                        arma::mat BsiKron = BCube.slice(i); // of dim O x hO
                        arma::mat Sigmasi = CovCube.slice(i); // of dim O x O
                        arma::colvec MeanJsi = BsiKron * LambdaOrigJNsi;
                        LambdaKrigJMat.col(i) = rmvnormRcpp(1, MeanJsi, Sigmasi);
                    }
                    LambdaKrigMat.col(j) = arma::reshape(LambdaKrigJMat, NNewLoc * O, 1);
                }
            }
            LambdaKrig.col(s) = arma::reshape(LambdaKrigMat, NNewLoc * O * K, 1);
            if (Verbose) {
                Rcpp::Rcout.precision(0);
                if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                    Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
            }
        }
    }    
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
    return LambdaKrig;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List AlphaKriggingFixedL(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {
    datobjPREDspat DatObj = ConvertDatObjPREDspatFixedL(DatObj_List);
    paraPREDspatFixedLCLlambda Para = ConvertParaPREDspatLambdaFixedLCL(Para_List);
    int T = DatObj.Nu;
    int K = DatObj.K;
    int M = DatObj.M;
    int O = DatObj.O;
    int L = DatObj.L;
    int IS = DatObj.IS;
    int SA = DatObj.SA;
    int h = DatObj.h;
    int alphaWsFiles = DatObj.alphaWsFiles;
    arma::Col<int> SeqL = DatObj.SeqL;
    int NNewLoc = DatObj.NNewLoc;
    arma::umat nnIndpred = DatObj.nnIndpred; // of dim NNewLoc x h
    arma::mat dist = DatObj.SpDist;
    arma::mat distOrigNew = DatObj.distOrigNew;
    arma::mat distNewNew = DatObj.distNewNew;
    // if IS = 1, then SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function
    int storeW = DatObj.storeW;
    arma::mat RhoMat = Para.Rho;
    arma::mat KappaMat = Para.Kappa;
    arma::mat ThetaMat = Para.Theta; // of dim NKeep x (K x L)
    arma::mat AlphaMat = Para.Alpha; // of dim NKeep x (M x O x K x (L - 1)) if alphaWsFiles = 0 and is not used (1 by 1 matrix of 1) if alphaWsFiles = 1
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Alpha: 0%.. ";
    }
    double rho;   
    arma::cube AlphaKrigCube(NNewLoc * O, L - 1, K); //rows of AlphaKrigCube first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on), same as that for AlphaOrig 
    arma::mat AlphaKrig(NKeep, K * (L - 1) * NNewLoc * O);
    //For each row (MCMC iteration) in AlphaKrig and WeightsKrig, the corresponding predicted entries for Alpha and Weights are ordered first by observation type, then spatially, 
    //by clustering group, and finally by factor. This is different from the ordering (for Alpha) for the outputs from the main function for the original M location points.
    arma::cube AlphaOrig(M * O, L - 1, K);
    if (O == 1) {
        double kappa;
        for (arma::uword s = 0; s < NKeep; s++) {
            rho = arma::as_scalar(RhoMat.row(s));
            kappa = arma::as_scalar(KappaMat.row(s));
            if (alphaWsFiles == 0) {
                arma::cube AlphaOrigVecCube(M * (L - 1) * K, 1, 1);
                AlphaOrigVecCube(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = arma::vectorise(AlphaMat.row(s));
                AlphaOrig = arma::reshape(AlphaOrigVecCube, M, L - 1, K);
            }
            else { // if alphaWsFiles = 1
                arma::cube AlphaOrigFile(L - 1, M, K);
                AlphaOrigFile.load("fixedLalphaKeptIter" + std::__cxx11::to_string(s + 1) + ".txt", arma::arma_ascii);
                AlphaOrig = arma::reshape(AlphaOrigFile, M, L - 1, K);
            }
            if (IS == 0) {
                arma::colvec Mean(NNewLoc, arma::fill::zeros);
                arma::mat FrhoNewNew(NNewLoc, NNewLoc, arma::fill::eye);
                arma::mat Sigma = kappa * FrhoNewNew;
                arma::cube AlphaKrigMatCube(NNewLoc, (L - 1) * K, 1);
                AlphaKrigMatCube(arma::span::all, arma::span::all, arma::span(0, 0)) = rmvnormRcpp(K * (L - 1), Mean, Sigma);
                AlphaKrigCube = arma::reshape(AlphaKrigMatCube, NNewLoc, L - 1, K);
            }
            else if (SA == 0) {// in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                //Calculate covariance components
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                arma::mat FrhoPlus = arma::trans(FrhoOrigNew) * CholInv(FrhoOrigOrig);
                arma::mat Sigma = kappa * (FrhoNewNew - FrhoPlus * FrhoOrigNew);
                for (arma::uword j = 0; j < K; j++) {
                    arma::mat AlphaOrigJ = AlphaOrig.slice(j); // of dim M x (L - 1) 
                    arma::mat AlphaKrigMatJ(NNewLoc, L - 1);
                    for (arma::uword l = 0; l < L - 1; l++) {
                        arma::colvec AlphaOrigJl = AlphaOrigJ.col(l);
                        arma::colvec MeanJl = FrhoPlus * AlphaOrigJl;
                        AlphaKrigMatJ.col(l) = rmvnormRcpp(1, MeanJl, Sigma);
                    }
                    AlphaKrigCube.slice(j) = AlphaKrigMatJ;
                }
            }
            else { //when include.space = spatApprox = TRUE; in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::vec sigma2vec(NNewLoc, arma::fill::ones);
                arma::mat Bmat(NNewLoc, h);
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    //arma::uvec nni = nnIndpred(arma::span(i, i), arma::span::all).t();
                    arma::uvec nni = arma::trans(nnIndpred.row(i));
                    arma::colvec FrhoOrigsi = FrhoOrigNew.col(i);
                    arma::colvec Frhosinsi = FrhoOrigsi(nni);//colvec
                    arma::mat Bsi = arma::trans(Frhosinsi) * CholInv(FrhoOrigOrig(nni, nni));
                    Bmat.row(i) = Bsi;
                    double Fsi = arma::as_scalar(1 - Bsi * Frhosinsi);
                    if (Fsi <= 0) Fsi = 0.00001;
                    sigma2vec(i) = Fsi * kappa;
                }      
                for (arma::uword j = 0; j < K; j++) {
                    arma::mat AlphaOrigJ = AlphaOrig.slice(j); // of dim M x (L - 1) 
                    arma::mat AlphaKrigMatJ(NNewLoc, L - 1);
                    for (arma::uword l = 0; l < L - 1; l++) {
                        arma::colvec AlphaOrigJl = AlphaOrigJ.col(l); // of length M
                        arma::colvec AlphaKrigJl(NNewLoc);
                        for (arma::uword i = 0; i < NNewLoc; i++) {
                            //arma::uvec nni = nnIndpred(arma::span(i, i), arma::span::all).t();
                            arma::uvec nni = arma::trans(nnIndpred.row(i));
                            arma::colvec AlphaOrigJlNsi = AlphaOrigJl(nni); // of length h
                            arma::mat Bsi = Bmat.row(i); // of dim 1 x h
                            double sdsi = sqrt(sigma2vec(i));
                            double meanJlsi = arma::as_scalar(Bsi * AlphaOrigJlNsi);
                            AlphaKrigJl(i) = arma::as_scalar(rnormRcpp(1, meanJlsi, sdsi));
                        }
                        AlphaKrigMatJ.col(l) = AlphaKrigJl;
                    }
                    AlphaKrigCube.slice(j) = AlphaKrigMatJ;
                }
            }
            AlphaKrig.row(s) = arma::trans(arma::vectorise(AlphaKrigCube));
            if (Verbose) {
                Rcpp::Rcout.precision(0);
                if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                    Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
            }
        }
    }
    else { // if O > 1
        arma::mat EyeO(O, O, arma::fill::eye);
        arma::rowvec KappaRow(O * (O + 1) / 2);
        arma::mat Kappa(O, O);
        int counter;
        for (arma::uword s = 0; s < NKeep; s++) {
            rho = arma::as_scalar(RhoMat.row(s));
            KappaRow = KappaMat.row(s);
            counter = 0;
            for (arma::uword i = 0; i < O; i++) {
                for (arma::uword j = 0; j <= i; j++) {
                    Kappa(i, j) = KappaRow(counter);
                    counter++;
                }
            }
            Kappa = arma::symmatl(Kappa);
            if (alphaWsFiles == 0) {
                arma::cube AlphaOrigVecCube(M * O * (L - 1) * K, 1, 1);
                AlphaOrigVecCube(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = arma::vectorise(AlphaMat.row(s));
                AlphaOrig = arma::reshape(AlphaOrigVecCube, M * O, L - 1, K);
            }
            else { // if alphaWsFiles = 1
                arma::cube AlphaOrigFile(L - 1, M * O, K);
                AlphaOrigFile.load("fixedLalphaKeptIter" + std::__cxx11::to_string(s + 1) + ".txt", arma::arma_ascii);
                AlphaOrig = arma::reshape(AlphaOrigFile, M * O, L - 1, K);
            }           
            if (IS == 0) {
                arma::colvec Mean(NNewLoc * O, arma::fill::zeros);
                arma::mat FrhoNewNew(NNewLoc, NNewLoc, arma::fill::eye);
                arma::mat Sigma = arma::kron(FrhoNewNew, Kappa);
                arma::cube AlphaKrigMatCube(NNewLoc * O, (L - 1) * K, 1);
                AlphaKrigMatCube(arma::span::all, arma::span::all, arma::span(0, 0)) = rmvnormRcpp(K * (L - 1), Mean, Sigma);
                AlphaKrigCube = arma::reshape(AlphaKrigMatCube, NNewLoc * O, L - 1, K);
            }
            else if (SA == 0) {// in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                //Calculate covariance components
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                arma::mat FrhoPlus = arma::trans(FrhoOrigNew) * CholInv(FrhoOrigOrig);
                arma::mat Sigma = arma::kron((FrhoNewNew - FrhoPlus * FrhoOrigNew), Kappa);
                arma::mat FrhoPlusKron = arma::kron(FrhoPlus, EyeO);
                for (arma::uword j = 0; j < K; j++) {
                    arma::mat AlphaOrigJ = AlphaOrig.slice(j); // of dim (M * O) x (L - 1) 
                    arma::mat AlphaKrigMatJ(NNewLoc * O, L - 1);
                    for (arma::uword l = 0; l < L - 1; l++) {
                        arma::colvec AlphaOrigJl = AlphaOrigJ.col(l);
                        arma::colvec MeanJl = FrhoPlusKron * AlphaOrigJl;
                        AlphaKrigMatJ.col(l) = rmvnormRcpp(1, MeanJl, Sigma);
                    }
                    AlphaKrigCube.slice(j) = AlphaKrigMatJ;
                }
            }
            else { //when include.space = spatApprox = TRUE; in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::cube CovCube(O, O, NNewLoc);
                arma::cube BCube(O, h * O, NNewLoc);
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    arma::uvec nni = arma::trans(nnIndpred.row(i));
                    arma::colvec FrhoOrigsi = FrhoOrigNew.col(i);
                    arma::mat Frhosinsi = FrhoOrigsi(nni);//colvec
                    arma::mat Bsi = arma::trans(Frhosinsi) * CholInv(FrhoOrigOrig(nni, nni));
                    double Fsi = arma::as_scalar(1 - Bsi * Frhosinsi);
                    if (Fsi <= 0) Fsi = 0.00001;
                    CovCube.slice(i) = Fsi * Kappa;
                    BCube.slice(i) = arma::kron(Bsi, EyeO);
                }
                for (arma::uword j = 0; j < K; j++) {
                    arma::mat AlphaOrigJ = AlphaOrig.slice(j); // of dim (M * O) x (L - 1) 
                    arma::mat AlphaKrigMatJ(NNewLoc * O, L - 1);
                    for (arma::uword l = 0; l < L - 1; l++) {
                        arma::colvec AlphaOrigJl = AlphaOrigJ.col(l);
                        arma::mat AlphaOrigJlMat = arma::reshape(AlphaOrigJl, O, M);
                        arma::mat AlphaKrigJlMat(O, NNewLoc);
                        for (arma::uword i = 0; i < NNewLoc; i++) {
                            arma::uvec nni = arma::trans(nnIndpred.row(i));
                            arma::mat AlphaOrigJlNsiMat = AlphaOrigJlMat.cols(nni); // of dim O x h
                            arma::colvec AlphaOrigJlNsi = arma::reshape(AlphaOrigJlNsiMat, O * h, 1);
                            arma::mat BsiKron = BCube.slice(i); // of dim O x hO
                            arma::mat Sigmasi = CovCube.slice(i); // of dim O x O
                            arma::colvec MeanJlsi = BsiKron * AlphaOrigJlNsi;
                            AlphaKrigJlMat.col(i) = rmvnormRcpp(1, MeanJlsi, Sigmasi);
                        }
                        AlphaKrigMatJ.col(l) = arma::reshape(AlphaKrigJlMat, NNewLoc * O, 1);
                    }
                    AlphaKrigCube.slice(j) = AlphaKrigMatJ;
                }
            }
            AlphaKrig.row(s) = arma::trans(arma::vectorise(AlphaKrigCube));
            if (Verbose) {
                Rcpp::Rcout.precision(0);
                if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                    Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
            }
        }
    }       
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

    // calculate weights and sample xi; then calculate the predicted factor loadings matrix Lambda
    arma::mat WeightsKrig(NKeep, K * L * NNewLoc * O, arma::fill::zeros);
    //For each row (MCMC iteration) in AlphaKrig and WeightsKrig, the corresponding predicted entries for Alpha and Weights are ordered first by observation type, then spatially, 
    //by clustering group, and finally by factor. This is different from the ordering (for Alpha) for the outputs from the main function for the original M location points.
    arma::mat LambdaKrig(NNewLoc * O * K, NKeep);
    //For each kept MCMC iteration, the corresponding predicted entries for the factor loadings matrix are ordered first by observation type, then spatially, and finally by factor. 
    //(the first (NNewLoc x O) entries correspond to factor 1, the next (NNewLoc x O) entries correspond to factor 2 and so on;
    //the first O rows correspond to the first new location point for prediction, the next O rows correspond to the second new location point for prediction and so on)                
    arma::cube WeightsKrigCube(NNewLoc * O, L, K, arma::fill::zeros); //rows of WeightsKrigCube first ordered by observation type and then spatially 
    //(the first O correspond to new loc 1, the next O correspond to new loc 2 and so on)
    arma::mat LambdaKrigMat(NNewLoc * O, K);//rows of LambdaKrigMat first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on), different as that for LambdaOrig
    arma::mat WeightsKrigMatJ(NNewLoc * O, L); //rows of WeightsKrigMatJ first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on)
    arma::colvec LambdaKrigJ(NNewLoc * O); //first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on), different as that for LambdaOrigJ 
    arma::mat Theta(L, K);
    arma::colvec ThetaJ(L);
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        Rcpp::Rcout << std::fixed << "Krigging Lambda: 0%.. ";
    }
    for (arma::uword s = 0; s < NKeep; s++) {
        Theta = arma::reshape(ThetaMat.row(s), L, K);
        arma::mat AlphaKrigS = arma::reshape(AlphaKrig.row(s), NNewLoc * O * (L - 1), K);
        for (arma::uword j = 0; j < K; j++) {
            ThetaJ = Theta.col(j); // of length L
            arma::mat AlphaKrigMatJ = arma::reshape(AlphaKrigS.col(j), NNewLoc * O, L - 1);
            for (arma::uword i = 0; i < NNewLoc; i++) {
                for (arma::uword o = 0; o < O; o++) {
                    arma::uword index = o + O * i;
                    double temp = 1;
                    for (arma::uword l = 0; l < L; l++) {
                        if (l == 0) { WeightsKrigMatJ(index, l) = pnormRcpp(AlphaKrigMatJ(index, l)); }
                        else if (l < L - 1) {
                            temp *= UpperpnormRcpp(AlphaKrigMatJ(index, l - 1));
                            WeightsKrigMatJ(index, l) = pnormRcpp(AlphaKrigMatJ(index, l)) * temp;
                        }
                        else {
                            WeightsKrigMatJ(index, l) = temp * UpperpnormRcpp(AlphaKrigMatJ(index, l - 1));
                        }
                    }
                    arma::vec ProbsJIO = WeightsKrigMatJ.row(index).t();
                    arma::uword XiJIO = arma::as_scalar(sampleRcpp(SeqL, 1, true, ProbsJIO));
                    LambdaKrigJ(index) = ThetaJ(XiJIO);
                }
            }
            WeightsKrigCube.slice(j) = WeightsKrigMatJ;
            LambdaKrigMat.col(j) = LambdaKrigJ;
        }
        WeightsKrig.row(s) = arma::trans(arma::vectorise(WeightsKrigCube));
        LambdaKrig.col(s) = arma::reshape(LambdaKrigMat, NNewLoc * O * K, 1);
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

    return Rcpp::List::create(Rcpp::Named("alpha") = AlphaKrig,
                              Rcpp::Named("weights") = WeightsKrig,
                              Rcpp::Named("lambda") = LambdaKrig);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List AlphaKriggingVaryLj(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {
    datobjPREDspat DatObj = ConvertDatObjPREDspatVaryLj(DatObj_List);
    int K = DatObj.K;
    paraPREDspatVaryLjLambda Para = ConvertParaPREDspatLambdaVaryLj(Para_List, K, NKeep);
    int T = DatObj.Nu;
    int M = DatObj.M;
    int O = DatObj.O;
    int IS = DatObj.IS;
    int SA = DatObj.SA;
    int h = DatObj.h;   
    int NNewLoc = DatObj.NNewLoc;
    arma::umat nnIndpred = DatObj.nnIndpred; // of dim NNewLoc x h
    arma::mat dist = DatObj.SpDist;
    arma::mat distOrigNew = DatObj.distOrigNew;
    arma::mat distNewNew = DatObj.distNewNew;
    // if IS = 1, then SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function
    arma::mat LjVecMat = Para.LjVec;
    arma::mat RhoMat = Para.Rho;
    arma::mat KappaMat = Para.Kappa;
    arma::field<arma::vec> ThetaField = Para.Theta; // the field is of dim K x NKeep, where each element vector is of length Lj (for that factor j at that MCMC iteration)
    arma::field<arma::mat> AlphaOrigField = Para.Alpha; // the field is of dim K x NKeep, where each element matrix is of dim (Lj - 1) x (M x O)
    //If Lj = 1 for a certain j at a particular MCMC iteration, then that corresponding matrix in alpha is of dimension 1 x (M x O) with entries all set to +Inf.
    //In each element matrix for each clustering group, the corresponding entries are ordered first by observation type, then spatially.
    arma::vec VerboseSeq;
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        VerboseSeq = floor(VerboseSeq);
        Rcpp::Rcout << std::fixed << "Krigging Alpha: 0%.. ";
    }
    double rho;    
    arma::field<arma::mat> AlphaKrigField(K, NKeep); //each element matrix is of dim (NNewLoc x O) x (Lj - 1) 
    //If Lj = 1 for a certain j at a particular MCMC iteration, then that corresponding matrix in alpha is of dimension (NNewLoc x O) x 1 with entries all set to +Inf.
    //In each element matrix for each clustering group, the corresponding entries are ordered first by observation type, then spatially.   
    if (O == 1) {
        double kappa;
        for (arma::uword s = 0; s < NKeep; s++) {
            rho = arma::as_scalar(RhoMat.row(s));
            kappa = arma::as_scalar(KappaMat.row(s));
            arma::colvec LjVec = LjVecMat.row(s).t();
            if (IS == 0) {
                arma::colvec Mean(NNewLoc, arma::fill::zeros);
                arma::mat FrhoNewNew(NNewLoc, NNewLoc, arma::fill::eye);
                arma::mat Sigma = kappa * FrhoNewNew;
                for (arma::uword j = 0; j < K; j++) {
                    int Lj = LjVec(j);
                    if (Lj == 1) {
                        arma::mat AlphaSJ(NNewLoc, 1);
                        AlphaSJ.fill(arma::datum::inf);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                    else { // Lj > 1
                        arma::mat AlphaSJ = rmvnormRcpp(Lj - 1, Mean, Sigma);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                }
            }
            else if (SA == 0) {// in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                arma::mat FrhoPlus = arma::trans(FrhoOrigNew) * CholInv(FrhoOrigOrig);
                arma::mat Sigma = kappa * (FrhoNewNew - FrhoPlus * FrhoOrigNew);
                for (arma::uword j = 0; j < K; j++) {
                    int Lj = LjVec(j);
                    if (Lj == 1) {
                        arma::mat AlphaSJ(NNewLoc, 1);
                        AlphaSJ.fill(arma::datum::inf);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                    else {
                        arma::mat AlphaOrigSJ = arma::trans(AlphaOrigField(j, s)); // of dim M x (Lj - 1) 
                        arma::mat AlphaKrigSJ(NNewLoc, Lj - 1);
                        for (arma::uword l = 0; l < Lj - 1; l++) {
                            arma::colvec AlphaOrigSJl = AlphaOrigSJ.col(l);
                            arma::colvec MeanSJl = FrhoPlus * AlphaOrigSJl;
                            AlphaKrigSJ.col(l) = rmvnormRcpp(1, MeanSJl, Sigma);
                        }
                        AlphaKrigField(j, s) = AlphaKrigSJ;
                    }
                }
            }
            else { //when include.space = spatApprox = TRUE; in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::vec sigma2vec(NNewLoc, arma::fill::ones);
                arma::mat Bmat(NNewLoc, h);
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    arma::uvec nni = arma::trans(nnIndpred.row(i));
                    arma::colvec FrhoOrigsi = FrhoOrigNew.col(i);
                    arma::colvec Frhosinsi = FrhoOrigsi(nni);//colvec
                    arma::mat Bsi = arma::trans(Frhosinsi) * CholInv(FrhoOrigOrig(nni, nni));
                    Bmat.row(i) = Bsi;
                    double Fsi = arma::as_scalar(1 - Bsi * Frhosinsi);
                    if (Fsi <= 0) {
                        //Rcpp::Rcout.precision(2);
                        //Rcpp::Rcout << "iter" << s << "newLoc" << i << "FsiValue" << Fsi << "\n";
                        Fsi = 0.00001;                       
                    }
                    sigma2vec(i) = Fsi * kappa;
                }
                for (arma::uword j = 0; j < K; j++) {
                    int Lj = LjVec(j);
                    if (Lj == 1) {
                        arma::mat AlphaSJ(NNewLoc, 1);
                        AlphaSJ.fill(arma::datum::inf);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                    else {
                        arma::mat AlphaOrigSJ = arma::trans(AlphaOrigField(j, s)); // of dim M x (Lj - 1) 
                        arma::mat AlphaKrigSJ(NNewLoc, Lj - 1);
                        for (arma::uword l = 0; l < Lj - 1; l++) {
                            arma::colvec AlphaOrigJl = AlphaOrigSJ.col(l); // of length M
                            arma::colvec AlphaKrigJl(NNewLoc);
                            for (arma::uword i = 0; i < NNewLoc; i++) {
                                arma::uvec nni = arma::trans(nnIndpred.row(i));
                                arma::colvec AlphaOrigJlNsi = AlphaOrigJl(nni); // of length h
                                arma::mat Bsi = Bmat.row(i); // of dim 1 x h
                                double sdsi = sqrt(sigma2vec(i));
                                double meanJlsi = arma::as_scalar(Bsi * AlphaOrigJlNsi);
                                AlphaKrigJl(i) = arma::as_scalar(rnormRcpp(1, meanJlsi, sdsi));
                            }
                            AlphaKrigSJ.col(l) = AlphaKrigJl;
                        }
                        AlphaKrigField(j, s) = AlphaKrigSJ;
                    }
                }                             
            }
            if (Verbose) {
                Rcpp::Rcout.precision(0);
                if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                    Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
            }
        }
    }
    else {
        arma::mat EyeO(O, O, arma::fill::eye);
        arma::rowvec KappaRow(O * (O + 1) / 2);
        arma::mat Kappa(O, O);
        int counter;
        for (arma::uword s = 0; s < NKeep; s++) {
            rho = arma::as_scalar(RhoMat.row(s));
            KappaRow = KappaMat.row(s);
            counter = 0;
            for (arma::uword i = 0; i < O; i++) {
                for (arma::uword j = 0; j <= i; j++) {
                    Kappa(i, j) = KappaRow(counter);
                    counter++;
                }
            }
            Kappa = arma::symmatl(Kappa);
            arma::colvec LjVec = LjVecMat.row(s).t();
            if (IS == 0) {
                arma::colvec Mean(NNewLoc * O, arma::fill::zeros);
                arma::mat FrhoNewNew(NNewLoc, NNewLoc, arma::fill::eye);
                arma::mat Sigma = arma::kron(FrhoNewNew, Kappa);
                for (arma::uword j = 0; j < K; j++) {
                    int Lj = LjVec(j);
                    if (Lj == 1) {
                        arma::mat AlphaSJ(NNewLoc * O, 1);
                        AlphaSJ.fill(arma::datum::inf);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                    else { // Lj > 1
                        arma::mat AlphaSJ = rmvnormRcpp(Lj - 1, Mean, Sigma);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                }
            }
            else if (SA == 0) {// in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                //Calculate covariance components
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                arma::mat FrhoPlus = arma::trans(FrhoOrigNew) * CholInv(FrhoOrigOrig);
                arma::mat Sigma = arma::kron((FrhoNewNew - FrhoPlus * FrhoOrigNew), Kappa);
                arma::mat FrhoPlusKron = arma::kron(FrhoPlus, EyeO);
                for (arma::uword j = 0; j < K; j++) {
                    int Lj = LjVec(j);
                    if (Lj == 1) {
                        arma::mat AlphaSJ(NNewLoc * O, 1);
                        AlphaSJ.fill(arma::datum::inf);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                    else {
                        arma::mat AlphaOrigSJ = arma::trans(AlphaOrigField(j, s)); // of dim (M * O) x (Lj - 1) 
                        arma::mat AlphaKrigSJ(NNewLoc * O, Lj - 1);
                        for (arma::uword l = 0; l < Lj - 1; l++) {
                            arma::colvec AlphaOrigSJl = AlphaOrigSJ.col(l);
                            arma::colvec MeanSJl = FrhoPlusKron * AlphaOrigSJl;
                            AlphaKrigSJ.col(l) = rmvnormRcpp(1, MeanSJl, Sigma);
                        }
                        AlphaKrigField(j, s) = AlphaKrigSJ;
                    }
                }
            }
            else { //when include.space = spatApprox = TRUE; in this case SpCorInd must = 0 (spatial.structure must be "continuous") as ensured by the corresponding R function 
                arma::cube CovCube(O, O, NNewLoc);
                arma::cube BCube(O, h * O, NNewLoc);
                arma::mat FrhoOrigOrig = SpEXP(rho, dist);
                arma::mat FrhoOrigNew = SpEXP(rho, distOrigNew);
                arma::mat FrhoNewNew = SpEXP(rho, distNewNew);
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    arma::uvec nni = arma::trans(nnIndpred.row(i));
                    arma::colvec FrhoOrigsi = FrhoOrigNew.col(i);
                    arma::mat Frhosinsi = FrhoOrigsi(nni);//colvec
                    arma::mat Bsi = arma::trans(Frhosinsi) * CholInv(FrhoOrigOrig(nni, nni));
                    double Fsi = arma::as_scalar(1 - Bsi * Frhosinsi);
                    if (Fsi <= 0) Fsi = 0.00001;
                    CovCube.slice(i) = Fsi * Kappa;
                    BCube.slice(i) = arma::kron(Bsi, EyeO);
                }
                for (arma::uword j = 0; j < K; j++) {
                    int Lj = LjVec(j);
                    if (Lj == 1) {
                        arma::mat AlphaSJ(NNewLoc * O, 1);
                        AlphaSJ.fill(arma::datum::inf);
                        AlphaKrigField(j, s) = AlphaSJ;
                    }
                    else {
                        arma::mat AlphaOrigSJ = arma::trans(AlphaOrigField(j, s)); // of dim (M * O) x (Lj - 1) 
                        arma::mat AlphaKrigSJ(NNewLoc * O, Lj - 1);
                        for (arma::uword l = 0; l < Lj - 1; l++) {
                            arma::colvec AlphaOrigSJl = AlphaOrigSJ.col(l);
                            arma::mat AlphaOrigSJlMat = arma::reshape(AlphaOrigSJl, O, M);
                            arma::mat AlphaKrigSJlMat(O, NNewLoc);
                            for (arma::uword i = 0; i < NNewLoc; i++) {
                                arma::uvec nni = arma::trans(nnIndpred.row(i));
                                arma::mat AlphaOrigSJlNsiMat = AlphaOrigSJlMat.cols(nni); // of dim O x h
                                arma::colvec AlphaOrigSJlNsi = arma::reshape(AlphaOrigSJlNsiMat, O * h, 1);
                                arma::mat BsiKron = BCube.slice(i); // of dim O x hO
                                arma::mat Sigmasi = CovCube.slice(i); // of dim O x O
                                arma::colvec MeanSJlsi = BsiKron * AlphaOrigSJlNsi;
                                AlphaKrigSJlMat.col(i) = rmvnormRcpp(1, MeanSJlsi, Sigmasi);
                            }
                            AlphaKrigSJ.col(l) = arma::reshape(AlphaKrigSJlMat, NNewLoc * O, 1);
                        }
                        AlphaKrigField(j, s) = AlphaKrigSJ;
                    }
                }
            }
            if (Verbose) {
                Rcpp::Rcout.precision(0);
                if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                    Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
            }
        }
    }    
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

    // calculate weights and sample xi; then calculate the predicted factor loadings matrix Lambda
    arma::field<arma::mat> WeightsKrigField(K, NKeep); //each element matrix is of dim (NNewLoc x O) x Lj
    //In each element matrix for each clustering group, the corresponding entries are ordered first by observation type, then spatially.   
    arma::mat LambdaKrig(NNewLoc * O * K, NKeep);
    //For each kept MCMC iteration, the corresponding predicted entries for the factor loadings matrix are ordered first by observation type, then spatially, and finally by factor. 
    //(the first (NNewLoc x O) entries correspond to factor 1, the next (NNewLoc x O) entries correspond to factor 2 and so on;
    //the first O rows correspond to the first new location point for prediction, the next O rows correspond to the second new location point for prediction and so on)                
    arma::mat LambdaKrigMat(NNewLoc * O, K);//rows of LambdaKrigMat first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on), different as that for LambdaOrig    
    arma::colvec LambdaKrigSJ(NNewLoc * O, arma::fill::ones); //first ordered by observation type and then spatially (the first O correspond to new loc 1, 
    //the next O correspond to new loc 2 and so on), different as that for LambdaOrigJ 
    if (Verbose) {
        VerboseSeq = { 0.25, 0.50, 0.75 }; 	// new standard modern C++ brace initialization
        VerboseSeq *= NKeep;
        Rcpp::Rcout << std::fixed << "Krigging Lambda: 0%.. ";
    }
    for (arma::uword s = 0; s < NKeep; s++) {
        arma::colvec LjVec = LjVecMat.row(s).t();
        for (arma::uword j = 0; j < K; j++) {
            int Lj = LjVec(j);           
            arma::vec ThetaSJ = ThetaField(j, s); // of length Lj
            if (Lj == 1) {
                LambdaKrigSJ += arma::as_scalar(ThetaSJ);
                LambdaKrigMat.col(j) = LambdaKrigSJ;
                arma::mat WeightsKrigMatSJ(NNewLoc * O, Lj, arma::fill::ones);
                WeightsKrigField(j, s) = WeightsKrigMatSJ;
            }
            else {
                arma::mat AlphaKrigSJ = AlphaKrigField(j, s); // of dim (NNewLoc * O) x (Lj - 1)
                arma::mat WeightsKrigMatSJ(NNewLoc * O, Lj); //rows of WeightsKrigMatSJ first ordered by observation type and then spatially 
                //(the first O correspond to new loc 1, the next O correspond to new loc 2 and so on)
                arma::Col<int> SeqLjVec(Lj);
                for (arma::uword seq = 0; seq < Lj; seq++) SeqLjVec(seq) = seq;
                for (arma::uword i = 0; i < NNewLoc; i++) {
                    for (arma::uword o = 0; o < O; o++) {
                        arma::uword index = o + O * i;
                        double temp = 1;
                        for (arma::uword l = 0; l < Lj; l++) {
                            if (l == 0) { WeightsKrigMatSJ(index, l) = pnormRcpp(AlphaKrigSJ(index, l)); }
                            else if (l < Lj - 1) {
                                temp *= UpperpnormRcpp(AlphaKrigSJ(index, l - 1));
                                WeightsKrigMatSJ(index, l) = pnormRcpp(AlphaKrigSJ(index, l)) * temp;
                            }
                            else {
                                WeightsKrigMatSJ(index, l) = temp * UpperpnormRcpp(AlphaKrigSJ(index, l - 1));
                            }
                        }
                        arma::vec ProbsJIO = WeightsKrigMatSJ.row(index).t();
                        arma::uword XiJIO = arma::as_scalar(sampleRcpp(SeqLjVec, 1, true, ProbsJIO));
                        LambdaKrigSJ(index) = ThetaSJ(XiJIO);
                    }
                }
                WeightsKrigField(j, s) = WeightsKrigMatSJ;
                LambdaKrigMat.col(j) = LambdaKrigSJ;
            }            
        }
        LambdaKrig.col(s) = arma::reshape(LambdaKrigMat, NNewLoc * O * K, 1);
        if (Verbose) {
            Rcpp::Rcout.precision(0);
            if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
                Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
        }
    }
    if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

    return Rcpp::List::create(Rcpp::Named("alpha") = AlphaKrigField,
           Rcpp::Named("weights") = WeightsKrigField,
           Rcpp::Named("lambda") = LambdaKrig);
}