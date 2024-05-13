#define ARMA_DONT_PRINT_ERRORS // To suppress the cholesky warning
#include <iostream>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include "MCMC_bfaSpatTemp.h"

//Sample from a categorical distribution---------------------------------------------------
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob) {
  Rcpp::IntegerVector xIV = Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(x));
  Rcpp::NumericVector probNV = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(prob));
  Rcpp::IntegerVector ret = Rcpp::RcppArmadillo::sample(xIV, size, replace, probNV);
  return Rcpp::as<arma::vec>(ret);
}



//Log density of a multivariate normal-----------------------------------------------------
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti){
  arma::vec Z = arma::vectorise(arma::trans(Rooti) * (Y - Mu));
  arma::vec Pi2(1); Pi2(0) = 2.0 * M_PI;
  return (-(Y.size() / 2.0) * arma::log(Pi2) - 0.5 * (arma::trans(Z) * Z) + arma::sum(arma::log(arma::diagvec(Rooti))))[0];
}



//Log density of a normal distribution
double dlnorm(double x, double mu, double sigma2) {
    return -0.5 * log(2.0 * M_PI * sigma2) - 0.5 * ((x - mu) * (x - mu)) / sigma2;
}



//Log density of a normal distribution - vectorized
arma::colvec dlnormRcpp(arma::vec const& x, arma::vec const& mu, arma::vec const& sigma2) {
  return -0.5 * log(2 * M_PI * sigma2) - 0.5 * ((x - mu) % (x - mu)) % (1 / sigma2);
}



//Function for sampling random standard uniform variable-----------------------------------
double randuRcpp() {
  return R::runif(0, 1);
}



//Function for sampling random binomial variable-----------------------------------
double rbinomRcpp(double n, double p) {
  return R::rbinom(n, p);
}



//Log density of binomial-----------------------------------
double dlbinom(int x, int n, double pi) { 
  double temp = std::lgamma(n + 1) - std::lgamma(x + 1) - std::lgamma((n - x) + 1);
  temp += x * std::log(pi) + (n - x) * std::log(1 - pi);
  return temp;
}



//Function for sampling random chi square variable-----------------------------------------
double rchisqRcpp(double df) {
  return R::rchisq(df);
}



//Distribution function of a standard normal-----------------------------------------------
double pnormRcpp(double q) {
  return R::pnorm5(q, 0, 1, 1, 0);
}



//Log distribution function of a standard normal-----------------------------------------------
double lpnormRcpp(double q) {
  return R::pnorm5(q, 0, 1, 1, 1);
}



//Distribution function of a standard normal upper tail-----------------------------------------------
double UpperpnormRcpp(double q) {
  return R::pnorm5(q, 0, 1, 0, 0);
}



//Log distribution function of a standard normal upper tail-----------------------------------------------
double lUpperpnormRcpp(double q) {
  return R::pnorm5(q, 0, 1, 0, 1);
}



//Inverse distribution function of a standard normal---------------------------------------
double qnormRcpp(double p) {
    return R::qnorm5(p, 0, 1, 1, 0);
}

//Inverse distribution function of a standard normal upper tail---------------------------------------
double UpperqnormRcpp(double p) {
    return R::qnorm5(p, 0, 1, 0, 0);
}



//Function for sampling from a standard normal distribution--------------------------------
arma::vec rnormSNRcpp(int n) {
  arma::vec out(n);
  for (int i = 0; i < n; i++) out(i) = R::rnorm(0, 1);
  return out;
}



//Sample from a standard normal distribution-----------------------------------------------
arma::vec rnormRcpp(int n, double mean, double sd) {
  arma::vec muvec(1);
  muvec(0) = mean;
  return arma::repmat(muvec, n, 1) + rnormSNRcpp(n) * sd;
}



//Sample from a multivariate normal distribution-------------------------------------------
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma) {
  int ncols = sigma.n_cols;
  arma::vec Yvec(n * ncols);
  Yvec = rnormSNRcpp(n * ncols);
  arma::mat Y = arma::reshape(Yvec, n, ncols);
  arma::mat cholSigma;
  try {
      cholSigma = arma::chol(sigma);
  } catch(...) {
      cholSigma = getCholRobust(sigma);
  }
  return arma::trans(arma::repmat(mean, 1, n).t() + Y * cholSigma);
}



//Sample from a standard normal distribution vectorize-----------------------------------------------
arma::vec rnormVecRcpp(arma::vec const& mean, arma::vec const& sd) {
  return mean + rnormSNRcpp(mean.size()) % sd;
}



//Sample from an inverse gamma distribution-------------------------------------------
double rigammaRcpp(double Alpha, double Theta) {
  return 1 / R::rgamma(Alpha, 1 / Theta);
}



//Sample from a gamma distribution-------------------------------------------
double rgammaRcpp(double Alpha, double Theta) {
  return R::rgamma(Alpha, 1 / Theta);
}



// Sample from a normal distribution truncated below or above by a real number x. -----------------------
// This first version uses inverse sampling from a CDF. It has problems
// when the data gets too far from x due to precision issues with qnorm and pnorm.
// See "Simulation of truncated normal variables" C. Robert (2009). As such we also
// included truncated normal from the msm R package. This implementation uses an
// accept-reject method that is more efficient than the inverse sampling version
// for values far from x.
double rtnormRcpp(double mean, double sd, bool Above, double x) {

  //Declare Variables
  double RandU = randuRcpp();
  double cutoff = (x - mean) / sd;

  //Truncation Above by cutoff
  if (Above) return sd * qnormRcpp(RandU * pnormRcpp(cutoff)) + mean;// the value is from negative infinity to cutoff

  //Truncation Below by cutoff
  else return sd * qnormRcpp(RandU - pnormRcpp(cutoff) * (RandU - 1)) + mean;// the value is from cutoff to positive infinity 

}

double rtnormRcppMSM(double mean, double sd, double lower, double upper) {

  //Set truncated normal function
  // Rcpp::Environment msm("package:msm"); //this is equivalent to library(msm)
  Rcpp::Environment msm = Rcpp::Environment::namespace_env("msm"); //This is equivalent to PACKAGE::FUNCTION()
  Rcpp::Function rtnorm = msm["rtnorm"];

  //Evaluate pmvnorm
  SEXP rtnormSEXP;
  rtnormSEXP = rtnorm(Rcpp::Named("n", 1),
                      Rcpp::Named("mean", mean),
                      Rcpp::Named("sd", sd),
                      Rcpp::Named("lower", lower),
                      Rcpp::Named("upper", upper));

  //Convert output to double
  return Rcpp::as<double>(rtnormSEXP);

}


// Sample from a truncated multivariate normal distribution
arma::rowvec rtmvnormRcpp(arma::vec const& mean, arma::mat const& sigma, arma::vec const& lower, arma::vec const& upper) {
    // Rcpp::Environment tmvtnorm("package:tmvtnorm"); //this is equivalent to library(tmvtnorm)
    Rcpp::Environment tmvtnorm = Rcpp::Environment::namespace_env("tmvtnorm");
    Rcpp::Function rtmvnorm = tmvtnorm["rtmvnorm"];
    
    Rcpp::NumericVector Mean = Rcpp::NumericVector(mean.begin(), mean.end());
    Rcpp::NumericVector Lower = Rcpp::NumericVector(lower.begin(), lower.end());
    Rcpp::NumericVector Upper = Rcpp::NumericVector(upper.begin(), upper.end());
    SEXP rtmvnormSEXP;
    rtmvnormSEXP = rtmvnorm(Rcpp::Named("n", 1),
                            Rcpp::Named("mean", Mean),
                            Rcpp::Named("sigma", sigma),
                            Rcpp::Named("lower", Lower),
                            Rcpp::Named("upper", Upper),
                            Rcpp::Named("algorithm", "rejection"));

    return Rcpp::as<arma::rowvec>(rtmvnormSEXP);
}



//Sample from a Wishart distribution using the Bartlett decomposition----------------------
arma::mat rwishRcpp(double n, arma::mat const& V) {
  int p = V.n_rows;
  arma::mat L;
  try {
      L = arma::chol(V);
  } catch(...) {
      L = getCholRobust(V);
  }
  arma::mat A(p, p, arma::fill::zeros);
  arma::vec XiSqr(1);
  for (int i = 0; i < p; i++) {
    XiSqr(0) = rchisqRcpp(n - i);
    A(i, i) = arma::as_scalar(arma::sqrt(XiSqr)); 
  }
  if (p > 1) {
    arma::vec RandSN = rnormSNRcpp(p * (p - 1) / 2);
    int counter = 0;
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < i; j++){
        A(j, i) = RandSN(counter);
        counter++;
      }
    }
  }
  arma::mat AL = A * L;
  return arma::trans(AL) * AL;
}


//Sample from an Inverse-Wishart distribution----------------------------------------------
//arma::mat riwishRcpp(double n, arma::mat const& V) {
//  return arma::inv_sympd(rwishRcpp(n, arma::inv_sympd(V)));
//}


//Vectorized Polya-Gamma via the pgdraw package------------------------------------------------------------------
arma::vec pgRcpp(arma::vec const& b, arma::vec const& c) {
  
  //This is equivalent to library(pgdraw)
  Rcpp::Environment pgdraw_env = Rcpp::Environment::namespace_env("pgdraw"); //This is equivalent to PACKAGE::FUNCTION()
  Rcpp::Function pgdraw = pgdraw_env["pgdraw"];
  
  Rcpp::NumericVector B = Rcpp::NumericVector(b.begin(), b.end());
  Rcpp::NumericVector C = Rcpp::NumericVector(c.begin(), c.end());
  SEXP pgdrawSEXP = pgdraw(Rcpp::Named("b", B),
                           Rcpp::Named("c", C));
  
  //Convert output to Armadillo vector
  return Rcpp::as<arma::vec>(pgdrawSEXP);
  
}
