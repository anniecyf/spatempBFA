#include <iostream>
#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"


//Function to update the LjVec vector-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec UpdateLjVec(arma::mat const& U, arma::field<arma::mat> Weights, int K, int M, int O, arma::colvec LjVec) {
    arma::colvec LjVecNew(K);
    arma::uword Index;
    for (arma::uword j = 0; j < K; j++) {
        if (LjVec(j) == 1) { LjVecNew(j) = 1; }
        else {
            arma::colvec UJ = U.col(j);
            arma::mat WeightsJ = Weights(j);
            arma::colvec LjOI(M * O, arma::fill::ones);//impt!!
            for (arma::uword o = 0; o < O; o++) {
                for (arma::uword i = 0; i < M; i++) {
                    Index = i + M * o;                   
                    double oneMinusUJOI = 1 - UJ(Index);
                    arma::uvec foundCl = arma::find(arma::cumsum(WeightsJ.col(Index)) > oneMinusUJOI, 1);
                    if (size(foundCl, 0) == 1) LjOI(Index) += arma::as_scalar(foundCl);
                    else LjOI(Index) = LjVec(j); //size(foundCl, 0) == 0
                    //LjOI(Index) += arma::as_scalar(arma::find(arma::cumsum(WeightsJ.col(Index)) > oneMinusUJOI, 1));
                }
            }
            LjVecNew(j) = arma::max(LjOI);
        }        
    }
    return LjVecNew;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::field<arma::mat> UpdateAlpha(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec) {
    arma::field<arma::mat> AlphaNew(K);
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj == 1) {
            arma::mat AlphaJNew(1, M * O);
            AlphaJNew.fill(arma::datum::inf);
            AlphaNew(j) = AlphaJNew;
        }
        else {
            arma::mat AlphaJ = Alpha(j); // of dim (Lj - 1) x MO, where Lj is the previous MCMC iteration's Lj, bigger than or equal to the Lj in this current MCMC iteration
            arma::mat AlphaJNew = AlphaJ(arma::span(0, Lj - 2), arma::span::all);  // of dim (Lj - 1, M * O) the current updated Lj
            AlphaNew(j) = AlphaJNew;
        }
    }
    return AlphaNew;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::field<arma::mat> GetWeightsVaryLj(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec) {
    arma::field<arma::mat> Weights(K);
    arma::uword Index1;
    arma::uword Index2;
    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        if (Lj == 1) { Weights(j) = arma::mat(1, M * O, arma::fill::ones); }
        else {
            arma::mat WeightsJ(Lj, M * O, arma::fill::zeros);
            arma::mat AlphaJ = Alpha(j);
            for (arma::uword o = 0; o < O; o++) {
                for (arma::uword i = 0; i < M; i++) {
                    Index1 = i + M * o;
                    Index2 = o + O * i;
                    double temp = 1;
                    for (arma::uword l = 0; l < Lj; l++) {
                        if (l == 0) { WeightsJ(l, Index1) = pnormRcpp(AlphaJ(l, Index2)); }
                        else if (l < Lj - 1) {
                            temp *= UpperpnormRcpp(AlphaJ(l - 1, Index2));
                            WeightsJ(l, Index1) = pnormRcpp(AlphaJ(l, Index2)) * temp;
                        }
                        else {
                            WeightsJ(l, Index1) = temp * UpperpnormRcpp(AlphaJ(l - 1, Index2));
                        }
                    }
                }
            }
            Weights(j) = WeightsJ;
        }
    }
    return Weights;
}


// Function to initialize Theta, Alpha, Weights, logWeights for the VaryingLjs main function
// [[Rcpp::depends(RcppArmadillo)]]
std::pair<listParaVaryLj, paraVaryLj> InitializeListPara(datobjVaryLj DatObj, paraVaryLj Para) {
    int O = DatObj.O;
    int M = DatObj.M;
    int K = DatObj.K;
    int T = DatObj.Nu;
    arma::colvec LjVec = Para.LjVec;
    arma::umat Xi = Para.Xi;
    arma::colvec Tau = Para.Tau;
    arma::colvec Eta = Para.Eta;
    arma::colvec XBeta = Para.XBeta;
    arma::mat EyeT(T, T, arma::fill::eye);
    arma::field<arma::vec> Theta(K);
    arma::field<arma::mat> Alpha(K);

    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        Theta(j) = rnormRcpp(Lj, 0, std::sqrt(1 / Tau(j)));
        if (Lj == 1) {
            arma::mat AlphaJ(1, M * O);
            AlphaJ.fill(arma::datum::inf);
            Alpha(j) = AlphaJ;
        }
        else {
            arma::mat AlphaJ(Lj - 1, M * O, arma::fill::zeros);
            Alpha(j) = AlphaJ;
        }       
        //arma::mat Weights(Lj , M * O, arma::fill::zeros);
    }
    arma::field<arma::mat> Weights = GetWeightsVaryLj(Alpha, K, M, O, LjVec);

    arma::mat Lambda(M * O, K);
    for (arma::uword j = 0; j < K; j++) {
        arma::vec ThetaJ = Theta(j);
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {               
                Lambda(i + M * o, j) = ThetaJ(Xi(i + M * o, j));
            }
        }
    }

    listParaVaryLj ListPara{};
    ListPara.Theta = Theta;
    ListPara.Alpha = Alpha;
    ListPara.Weights = Weights;
    Para.Lambda = Lambda;
    Para.Mean = arma::kron(EyeT, Lambda) * Eta + XBeta;
    return std::pair<listParaVaryLj, paraVaryLj>(ListPara, Para);
}

std::pair<listParaVaryLj, VAR1paraVaryLj> InitializeListPara(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para) {
    int O = DatObj.O;
    int M = DatObj.M;
    int K = DatObj.K;
    int T = DatObj.Nu;
    arma::colvec LjVec = Para.LjVec;
    arma::umat Xi = Para.Xi;
    arma::colvec Tau = Para.Tau;
    arma::colvec Eta = Para.Eta;
    arma::colvec XBeta = Para.XBeta;
    arma::mat EyeT(T, T, arma::fill::eye);
    arma::field<arma::vec> Theta(K);
    arma::field<arma::mat> Alpha(K);

    for (arma::uword j = 0; j < K; j++) {
        int Lj = LjVec(j);
        Theta(j) = rnormRcpp(Lj, 0, std::sqrt(1 / Tau(j)));
        if (Lj == 1) {
            arma::mat AlphaJ(1, M * O);
            AlphaJ.fill(arma::datum::inf);
            Alpha(j) = AlphaJ;
        }
        else {
            arma::mat AlphaJ(Lj - 1, M * O, arma::fill::zeros);
            Alpha(j) = AlphaJ;
        }       
        //arma::mat Weights(Lj , M * O, arma::fill::zeros);
    }
    arma::field<arma::mat> Weights = GetWeightsVaryLj(Alpha, K, M, O, LjVec);

    arma::mat Lambda(M * O, K);
    for (arma::uword j = 0; j < K; j++) {
        arma::vec ThetaJ = Theta(j);
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {               
                Lambda(i + M * o, j) = ThetaJ(Xi(i + M * o, j));
            }
        }
    }

    listParaVaryLj ListPara{};
    ListPara.Theta = Theta;
    ListPara.Alpha = Alpha;
    ListPara.Weights = Weights;
    Para.Lambda = Lambda;
    Para.Mean = arma::kron(EyeT, Lambda) * Eta + XBeta;
    return std::pair<listParaVaryLj, VAR1paraVaryLj>(ListPara, Para);
}


//Functions to get weights, w_jl(s_i) and log weights, log(w_jl(s_i)) -------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube GetWeightsFixedL(arma::cube const& Alpha, int K, int M, int L, int O) {
  arma::cube Weights(L, M * O, K, arma::fill::zeros);
  arma::uword Index1;
  arma::uword Index2;
  for (arma::uword j = 0; j < K; j++) {
      for (arma::uword o = 0; o < O; o++) {
          for (arma::uword i = 0; i < M; i++) {
              Index1 = i + M * o;
              Index2 = o + O * i;
              double temp = 1;
              for (arma::uword l = 0; l < L; l++) {
                  if (l == 0) { Weights(l, Index1, j) = pnormRcpp(Alpha(l, Index2, j)); }
                  else if (l < L - 1) {
                      temp *= UpperpnormRcpp(Alpha(l - 1, Index2, j));
                      Weights(l, Index1, j) = pnormRcpp(Alpha(l, Index2, j)) * temp;
                  }
                  else {
                      Weights(l, Index1, j) = temp * UpperpnormRcpp(Alpha(l - 1, Index2, j));
                  }
              }
          }
      }
  }
  return Weights;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube GetLogWeightsFixedL(arma::cube const& Alpha, int K, int M, int L, int O) {
    arma::cube logWeights(L, M * O, K);
    logWeights.fill(-arma::datum::inf);
    arma::uword Index1;
    arma::uword Index2;
    for (arma::uword j = 0; j < K; j++) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index1 = i + M * o;
                Index2 = o + O * i;
                double temp = 0;
                for (arma::uword l = 0; l < L; l++) {
                    if (l == 0) { logWeights(l, Index1, j) = lpnormRcpp(Alpha(l, Index2, j)); }
                    else if (l < L - 1) {
                        temp += lUpperpnormRcpp(Alpha(l - 1, Index2, j));
                        logWeights(l, Index1, j) = lpnormRcpp(Alpha(l, Index2, j)) + temp;
                    }
                    else {
                        logWeights(l, Index1, j) = temp + lUpperpnormRcpp(Alpha(l - 1, Index2, j));
                    }

                }
            }
        }
    }
    return logWeights;
}

 
//Matrix inversion using cholesky decomposition for covariances-------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CholInv(arma::mat const& Cov) {
    arma::mat CovInv;
    if (Cov.is_sympd()) {
        try {
            CovInv = arma::inv_sympd(Cov);
        } catch(...) {
            arma::mat Z = getCloseNNDmat(Cov);
            int nrow = size(Cov, 0);
            arma::mat diagEye(nrow, nrow, arma::fill::eye);
            CovInv = arma::inv_sympd(Z + 0.0001 * diagEye);
        }
    } 
    else {
        try {
            CovInv = arma::inv(Cov);
        } catch(...) {
            arma::mat Z = getCloseNNDmat(Cov);
            int nrow = size(Cov, 0);
            arma::mat diagEye(nrow, nrow, arma::fill::eye);
            CovInv = arma::inv(Z + 0.0001 * diagEye);
        }
    }
    return CovInv;
}



// function used when scenario = 3 for the FixedL main function
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]] 
arma::field<arma::mat> whichJsCalc(arma::umat nnInd) {
    int M = size(nnInd, 0);
    int h = size(nnInd, 1);
    arma::field<arma::mat> whichJs(M, 1);
    int nnMaxNum;
    for (int i = 0; i < M; i++) {
        arma::mat iAsNN(2, 1, arma::fill::zeros);
        if (i < M - 1) {
            int iNumj = 0;
            for (int j = i + 1; j < M; j++) {
                if (j < h) nnMaxNum = j - 1;
                else nnMaxNum = h - 1;
                arma::uvec nnj = nnInd(arma::span(j, j), arma::span(0, nnMaxNum)).t();
                arma::uvec findVec = find(nnj == i);
                if (findVec.size() > 0) {
                    int ksi = findVec(0);
                    arma::vec ijInfo = { j, ksi };
                    if (iNumj == 0) iAsNN.col(0) = ijInfo;
                    else iAsNN.insert_cols(iNumj, ijInfo);
                    iNumj += 1;
                }
            }
        }
        whichJs(i, 0) = iAsNN;
    }
    return whichJs;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getCloseNNDmat(arma::mat matx){
	arma::mat Y = 0.5 * (matx + matx.t());
	int nrow = size(matx, 0);
	arma::vec vec0(nrow, arma::fill::zeros);
	arma::vec eigvalVec;
    arma::mat eigvecMat;
    arma::eig_sym(eigvalVec, eigvecMat, Y);
	arma::vec eigvalVecNN = arma::max(eigvalVec, vec0);
	arma::mat Z = eigvecMat * diagmat(eigvalVecNN) * eigvecMat.t();
	return Z;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getCholRobust(arma::mat matx){    
    arma::mat cholRobust;
    bool cholSuccess = false;
    while (!cholSuccess) {
        cholSuccess = arma::chol(cholRobust, matx); 
        if (cholSuccess == false){          
            arma::mat Z = getCloseNNDmat(matx);
            cholSuccess = arma::chol(cholRobust, Z);  
            int nrow = size(matx, 0);
            arma::mat diagEye(nrow, nrow, arma::fill::eye);
            while (!cholSuccess){
                Z += 0.0001 * diagEye;
                cholSuccess = arma::chol(cholRobust, Z); 
            }
        }
    } 
    return cholRobust;
}



//Function that checks numerical equality of two objects against a tolerance-----------------------------
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol) {
    return arma::approx_equal(lhs, rhs, "absdiff", tol);
}
