#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"

//Calculate inverse root of upper cholesky for log density of normal-----------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat GetRooti(arma::mat const& Cov) {
    arma::mat cholCov;
    try {
        cholCov = arma::chol(Cov);
    } catch(...) {
        cholCov = getCholRobust(Cov);
    }
    return arma::inv(arma::trimatu(cholCov));
}

// Only for the case when all time points are equi-spaced
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getRootiH(double Psi, int Nu, int TempCorInd, int seasonPeriod){
    double psi = Psi;//default AR(1), SAR(1)
    if (TempCorInd == 0 || TempCorInd == 3) {psi = exp(-Psi);}// for "exponential" or "sexponential" temporal.structure
    arma::mat rootiH(Nu, Nu, arma::fill::zeros);
    if (TempCorInd == 0 || TempCorInd == 1) {//no seasonality
        rootiH.diag() += (1 / sqrt(1 - std::pow(psi, 2)));
        rootiH(0, 0) = 1;
        rootiH.diag(1) -= psi / sqrt(1 - std::pow(psi, 2));
    }
    else if (TempCorInd == 2 || TempCorInd == 3) {//with seasonality
        rootiH.diag() += (1 / sqrt(1 - std::pow(psi, 2)));
        rootiH(arma::span(0, seasonPeriod - 1), arma::span(0, seasonPeriod - 1)).diag() *= sqrt(1 - std::pow(psi, 2));
        rootiH.diag(seasonPeriod) -= psi / sqrt(1 - std::pow(psi, 2));
    }
    
    return rootiH;
}

// Only for the case when all time points are equi-spaced 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getInvH(double Psi, int Nu, int TempCorInd, int seasonPeriod){
  double psi = Psi;//default AR(1), SAR(1)
  if (TempCorInd == 0 || TempCorInd == 3) { psi = exp(-Psi); }// for "exponential" or "sexponential" temporal.structure
  arma::mat invH(Nu, Nu, arma::fill::zeros);
  if (TempCorInd == 0 || TempCorInd == 1) {//no seasonality
      invH.diag() += 1 + std::pow(psi, 2);
      invH(0, 0) = 1;
      invH((Nu - 1), (Nu - 1)) = 1;
      invH.diag(1) -= psi;
      invH.diag(-1) -= psi;
      invH *= 1 / (1 - std::pow(psi, 2));
  }
  else if (TempCorInd == 2 || TempCorInd == 3) {//with seasonality
      invH.diag() += 1 + std::pow(psi, 2);
      invH(arma::span(0, seasonPeriod - 1), arma::span(0, seasonPeriod - 1)).diag() /= (1 + std::pow(psi, 2));
      invH(arma::span(Nu - seasonPeriod, Nu - 1), arma::span(Nu - seasonPeriod, Nu - 1)).diag() /= (1 + std::pow(psi, 2));
      invH.diag(seasonPeriod) -= psi;
      invH.diag(-seasonPeriod) -= psi;
      invH *= 1 / (1 - std::pow(psi, 2));
  }
  
  return invH;
}


//Function to calculate temporal correlation structure-------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat getH(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu, int seasonPeriod) {

  //Create ouput object
  arma::mat out(Nu, Nu);

  //exponential
  if (TempCorInd == 0) out = arma::exp(-Psi * TimeDist);

  //ar(1) continuous
  else if (TempCorInd == 1) {
      out = arma::eye(Nu, Nu);
      for (int j = 0; j < Nu; j++) {
          for (int k = 0; k < j; k++) {
              out(j, k) = std::pow(Psi, TimeDist(j, k));
          }
      }
      out = arma::symmatl(out);
  }

  //SAR(1)_seasonPeriod
  else if (TempCorInd == 2) {
      out = arma::eye(Nu, Nu);
      for (int j = 0; j < Nu; j++) {
          for (int k = 0; k < j; k++) {
              int power = TimeDist(j, k);
              if (power % seasonPeriod != 0) {
                  out(j, k) = 0;
              }
              else {
                  out(j, k) = std::pow(Psi, power / seasonPeriod);
              }
          }
      }
      out = arma::symmatl(out);
  }

  //Sexponential_seasonPeriod
  else if (TempCorInd == 3) {
      out = arma::eye(Nu, Nu);
      for (int j = 0; j < Nu; j++) {
          for (int k = 0; k < j; k++) {
              int power = TimeDist(j, k);
              if (power % seasonPeriod != 0) {
                  out(j, k) = 0;
              }
              else {
                  out(j, k) = exp(-Psi * (power / seasonPeriod));
              }
          }
      }
      out = arma::symmatl(out);
  }

  return out;
}



//Function to calculate spatial correlation structure-------------------------------------------------------------
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat SpEXP(double rho, arma::mat const& SpDist) {
    return arma::exp(-rho * SpDist);
}