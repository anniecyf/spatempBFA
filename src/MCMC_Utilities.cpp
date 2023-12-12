#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive) {

  //Set MCMC object
  int BarLength = McmcObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  0%..  ";
  }

}



//Function to pilot adapt tuning parameter--------------------------------------------------------------
double PilotAdaptFunc(double TuningParameter, double AcceptancePct) {

  //Adjust tuning parameter using scaling based on size of acceptance rate
  if (AcceptancePct >= 0.90) TuningParameter *= 1.3;
  if ( (AcceptancePct >= 0.75 ) & (AcceptancePct < 0.90 ) ) TuningParameter *= 1.2;
  if ( (AcceptancePct >= 0.45 ) & (AcceptancePct < 0.75 ) ) TuningParameter *= 1.1;
  if ( (AcceptancePct <= 0.25 ) & (AcceptancePct > 0.15 ) ) TuningParameter *= 0.9;
  if ( (AcceptancePct <= 0.15 ) & (AcceptancePct > 0.10 ) ) TuningParameter *= 0.8;
  if (AcceptancePct <= 0.10) TuningParameter *= 0.7;
  return TuningParameter;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj PilotAdaptation(metrobj MetrObj, datobjFixedL DatObj, mcmcobj McmcObj, int s) {
  
  //Set Metropolis objects
  double MetropPsi = MetrObj.MetropPsi;
  double AcceptancePsi = MetrObj.AcceptancePsi;

  //Get acceptance percentages
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;
  double PctPsi = AcceptancePsi / double(PilotAdaptDenominator);
  
  //Update Tuning Parameter
  MetropPsi = PilotAdaptFunc(MetropPsi, PctPsi);
  MetrObj.MetropPsi = MetropPsi;
  
  //Zero the acceptance counters
  AcceptancePsi = 0;
  MetrObj.AcceptancePsi = AcceptancePsi;
  
  //Update Rho
  if (DatObj.IS == 1) {
    double MetropRho = MetrObj.MetropRho;
    double AcceptanceRho = MetrObj.AcceptanceRho;
    double PctRho = AcceptanceRho / double(PilotAdaptDenominator);
    MetropRho = PilotAdaptFunc(MetropRho, PctRho);
    MetrObj.MetropRho = MetropRho;
    AcceptanceRho = 0;
    MetrObj.AcceptanceRho = AcceptanceRho;
  }  
  return MetrObj;

}

metrobj PilotAdaptation(metrobj MetrObj, datobjVaryLj DatObj, mcmcobj McmcObj, int s) {

    //Set Metropolis objects
    double MetropPsi = MetrObj.MetropPsi;
    double AcceptancePsi = MetrObj.AcceptancePsi;

    //Get acceptance percentages
    int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;
    double PctPsi = AcceptancePsi / double(PilotAdaptDenominator);

    //Update Tuning Parameter
    MetropPsi = PilotAdaptFunc(MetropPsi, PctPsi);
    MetrObj.MetropPsi = MetropPsi;

    //Zero the acceptance counters
    AcceptancePsi = 0;
    MetrObj.AcceptancePsi = AcceptancePsi;

    //Update Rho
    if (DatObj.IS == 1) {
        double MetropRho = MetrObj.MetropRho;
        double AcceptanceRho = MetrObj.AcceptanceRho;
        double PctRho = AcceptanceRho / double(PilotAdaptDenominator);
        MetropRho = PilotAdaptFunc(MetropRho, PctRho);
        MetrObj.MetropRho = MetropRho;
        AcceptanceRho = 0;
        MetrObj.AcceptanceRho = AcceptanceRho;
    }
    return MetrObj;

}


VAR1metrobj PilotAdaptation(VAR1metrobj MetrObj, VAR1datobjFixedL DatObj, mcmcobj McmcObj, int s) {
    
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;
  if (DatObj.IS == 1) {
    double MetropRho = MetrObj.MetropRho;
    double AcceptanceRho = MetrObj.AcceptanceRho;
    double PctRho = AcceptanceRho / double(PilotAdaptDenominator);
    MetropRho = PilotAdaptFunc(MetropRho, PctRho);
    MetrObj.MetropRho = MetropRho;
    AcceptanceRho = 0;
    MetrObj.AcceptanceRho = AcceptanceRho;
  }  
  return MetrObj;

}

VAR1metrobj PilotAdaptation(VAR1metrobj MetrObj, VAR1datobjVaryLj DatObj, mcmcobj McmcObj, int s) {

    int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;
    if (DatObj.IS == 1) {
        double MetropRho = MetrObj.MetropRho;
        double AcceptanceRho = MetrObj.AcceptanceRho;
        double PctRho = AcceptanceRho / double(PilotAdaptDenominator);
        MetropRho = PilotAdaptFunc(MetropRho, PctRho);
        MetrObj.MetropRho = MetropRho;
        AcceptanceRho = 0;
        MetrObj.AcceptanceRho = AcceptanceRho;
    }
    return MetrObj;

}




//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(metrobj MetrObj) {

  Rcpp::List Out;
  Out = Rcpp::List::create(Rcpp::Named("AcceptancePsi") = MetrObj.AcceptancePsi,
                           Rcpp::Named("MetropPsi") = MetrObj.MetropPsi,
                           Rcpp::Named("AcceptanceRho") = MetrObj.AcceptanceRho,
                           Rcpp::Named("MetropRho") = MetrObj.MetropRho);
  return Out;

}

Rcpp::List OutputMetrObj(VAR1metrobj MetrObj) {

  Rcpp::List Out;
  Out = Rcpp::List::create(Rcpp::Named("AcceptanceRho") = MetrObj.AcceptanceRho,
                           Rcpp::Named("MetropRho") = MetrObj.MetropRho);
  return Out;

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int s, mcmcobj McmcObj) {

  //Set MCMC object
  int NSims = McmcObj.NSims;
  int NBurn = McmcObj.NBurn;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  if (s < NSims + NBurn) Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%.. ";
  if (s == NSims + NBurn) Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%!";

}



//Functions for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, spatpara1 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    
    return col;
}

arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    int L = DatObj.L;
    int spatPred = DatObj.spatPred;
    int storeW = DatObj.storeW;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = ParaCL.Delta;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    arma::umat Xi = ParaCL.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::cube Alpha = ParaCL.Alpha;
    arma::cube Weights = ParaCL.Weights;
    arma::mat Theta = ParaCL.Theta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    if (spatPred == 1) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                for (arma::uword j = 0; j < K; j++) {
                    col(counter) = Xi(Index, j);
                    counter++;
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = Delta(j);
            counter++;
        }
        if (DatObj.alphaWsFiles == 0) {
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                for (arma::uword l = 0; l < (L - 1); l++) {
                    for (arma::uword i = 0; i < M; i++) {
                        for (arma::uword o = 0; o < O; o++) {
                            Index = o + O * i;
                            col(counter) = AlphaJ(l, Index);
                            counter++;
                        }
                    }
                }
            }
        }       
        for (arma::uword j = 0; j < K; j++) {
            for (arma::uword l = 0; l < L; l++) {
                col(counter) = Theta(l, j);
                counter++;
            }
        }
    }
    if (storeW == 1 && DatObj.alphaWsFiles == 0) {
        for (arma::uword j = 0; j < K; j++) {
            arma::mat WeightsJ = Weights.slice(j);
            for (arma::uword l = 0; l < L; l++) {
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index = i + M * o;
                        col(counter) = WeightsJ(l, Index);
                        counter++;
                    }
                }
            }
        }
    }

    return col;
}

arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, spatpara2 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }

    return col;
}

arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    int L = DatObj.L;
    int spatPred = DatObj.spatPred;
    int storeW = DatObj.storeW;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = ParaCL.Delta;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    arma::umat Xi = ParaCL.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::cube Alpha = ParaCL.Alpha;
    arma::cube Weights = ParaCL.Weights;
    arma::mat Theta = ParaCL.Theta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    if (spatPred == 1) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                for (arma::uword j = 0; j < K; j++) {
                    col(counter) = Xi(Index, j);
                    counter++;
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = Delta(j);
            counter++;
        }
        if (DatObj.alphaWsFiles == 0) {
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                for (arma::uword l = 0; l < (L - 1); l++) {
                    for (arma::uword i = 0; i < M; i++) {
                        for (arma::uword o = 0; o < O; o++) {
                            Index = o + O * i;
                            col(counter) = AlphaJ(l, Index);
                            counter++;
                        }
                    }
                }
            }
        }       
        for (arma::uword j = 0; j < K; j++) {
            for (arma::uword l = 0; l < L; l++) {
                col(counter) = Theta(l, j);
                counter++;
            }
        }
    }
    if (storeW == 1 && DatObj.alphaWsFiles == 0) {
        for (arma::uword j = 0; j < K; j++) {
            arma::mat WeightsJ = Weights.slice(j);
            for (arma::uword l = 0; l < L; l++) {
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index = i + M * o;
                        col(counter) = WeightsJ(l, Index);
                        counter++;
                    }
                }
            }
        }
    }

    return col;
}

arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, spatpara3 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }

    return col;
}

arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    int L = DatObj.L;
    int spatPred = DatObj.spatPred;
    int storeW = DatObj.storeW;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = ParaCL.Delta;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    arma::umat Xi = ParaCL.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::cube Alpha = ParaCL.Alpha;
    arma::cube Weights = ParaCL.Weights;
    arma::mat Theta = ParaCL.Theta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    if (spatPred == 1) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                for (arma::uword j = 0; j < K; j++) {
                    col(counter) = Xi(Index, j);
                    counter++;
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = Delta(j);
            counter++;
        }
        if (DatObj.alphaWsFiles == 0) {
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                for (arma::uword l = 0; l < (L - 1); l++) {
                    for (arma::uword i = 0; i < M; i++) {
                        for (arma::uword o = 0; o < O; o++) {
                            Index = o + O * i;
                            col(counter) = AlphaJ(l, Index);
                            counter++;
                        }
                    }
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            for (arma::uword l = 0; l < L; l++) {
                col(counter) = Theta(l, j);
                counter++;
            }
        }
    }
    if (storeW == 1 && DatObj.alphaWsFiles == 0) {
        for (arma::uword j = 0; j < K; j++) {
            arma::mat WeightsJ = Weights.slice(j);
            for (arma::uword l = 0; l < L; l++) {
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index = i + M * o;
                        col(counter) = WeightsJ(l, Index);
                        counter++;
                    }
                }
            }
        }
    }

    return col;
}


arma::colvec StoreSamples(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = Para.Delta;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    arma::umat Xi = Para.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::colvec LjVec = Para.LjVec;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Xi(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = Delta(j);
        counter++;
    }    
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = LjVec(j);
        counter++;
    }
    return col;
}

arma::colvec StoreSamples(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = Para.Delta;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    arma::umat Xi = Para.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::colvec LjVec = Para.LjVec;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Xi(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = Delta(j);
        counter++;
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = LjVec(j);
        counter++;
    }
    return col;
}

arma::colvec StoreSamples(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = Para.Delta;
    arma::mat Upsilon = Para.Upsilon;
    double Psi = Para.Psi;
    arma::umat Xi = Para.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::colvec LjVec = Para.LjVec;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    col(counter) = Psi;
    counter++;
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Xi(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = Delta(j);
        counter++;
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = LjVec(j);
        counter++;
    }
    return col;
}



arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara1 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    
    return col;
}

arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    int L = DatObj.L;
    int spatPred = DatObj.spatPred;
    int storeW = DatObj.storeW;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = ParaCL.Delta;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    arma::umat Xi = ParaCL.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::cube Alpha = ParaCL.Alpha;
    arma::cube Weights = ParaCL.Weights;
    arma::mat Theta = ParaCL.Theta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    if (spatPred == 1) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                for (arma::uword j = 0; j < K; j++) {
                    col(counter) = Xi(Index, j);
                    counter++;
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = Delta(j);
            counter++;
        }
        if (DatObj.alphaWsFiles == 0) {
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                for (arma::uword l = 0; l < (L - 1); l++) {
                    for (arma::uword i = 0; i < M; i++) {
                        for (arma::uword o = 0; o < O; o++) {
                            Index = o + O * i;
                            col(counter) = AlphaJ(l, Index);
                            counter++;
                        }
                    }
                }
            }
        }       
        for (arma::uword j = 0; j < K; j++) {
            for (arma::uword l = 0; l < L; l++) {
                col(counter) = Theta(l, j);
                counter++;
            }
        }
    }
    if (storeW == 1 && DatObj.alphaWsFiles == 0) {
        for (arma::uword j = 0; j < K; j++) {
            arma::mat WeightsJ = Weights.slice(j);
            for (arma::uword l = 0; l < L; l++) {
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index = i + M * o;
                        col(counter) = WeightsJ(l, Index);
                        counter++;
                    }
                }
            }
        }
    }

    return col;
}

arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara2 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }

    return col;
}

arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    int L = DatObj.L;
    int spatPred = DatObj.spatPred;
    int storeW = DatObj.storeW;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = ParaCL.Delta;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    arma::umat Xi = ParaCL.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::cube Alpha = ParaCL.Alpha;
    arma::cube Weights = ParaCL.Weights;
    arma::mat Theta = ParaCL.Theta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    if (spatPred == 1) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                for (arma::uword j = 0; j < K; j++) {
                    col(counter) = Xi(Index, j);
                    counter++;
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = Delta(j);
            counter++;
        }
        if (DatObj.alphaWsFiles == 0) {
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                for (arma::uword l = 0; l < (L - 1); l++) {
                    for (arma::uword i = 0; i < M; i++) {
                        for (arma::uword o = 0; o < O; o++) {
                            Index = o + O * i;
                            col(counter) = AlphaJ(l, Index);
                            counter++;
                        }
                    }
                }
            }
        }       
        for (arma::uword j = 0; j < K; j++) {
            for (arma::uword l = 0; l < L; l++) {
                col(counter) = Theta(l, j);
                counter++;
            }
        }
    }
    if (storeW == 1 && DatObj.alphaWsFiles == 0) {
        for (arma::uword j = 0; j < K; j++) {
            arma::mat WeightsJ = Weights.slice(j);
            for (arma::uword l = 0; l < L; l++) {
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index = i + M * o;
                        col(counter) = WeightsJ(l, Index);
                        counter++;
                    }
                }
            }
        }
    }

    return col;
}

arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara3 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }

    return col;
}

arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    int L = DatObj.L;
    int spatPred = DatObj.spatPred;
    int storeW = DatObj.storeW;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = ParaCL.Delta;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    arma::umat Xi = ParaCL.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::cube Alpha = ParaCL.Alpha;
    arma::cube Weights = ParaCL.Weights;
    arma::mat Theta = ParaCL.Theta;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    if (spatPred == 1) {
        for (arma::uword o = 0; o < O; o++) {
            for (arma::uword i = 0; i < M; i++) {
                Index = i + M * o;
                for (arma::uword j = 0; j < K; j++) {
                    col(counter) = Xi(Index, j);
                    counter++;
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = Delta(j);
            counter++;
        }
        if (DatObj.alphaWsFiles == 0) {
            for (arma::uword j = 0; j < K; j++) {
                arma::mat AlphaJ = Alpha.slice(j);
                for (arma::uword l = 0; l < (L - 1); l++) {
                    for (arma::uword i = 0; i < M; i++) {
                        for (arma::uword o = 0; o < O; o++) {
                            Index = o + O * i;
                            col(counter) = AlphaJ(l, Index);
                            counter++;
                        }
                    }
                }
            }
        }
        for (arma::uword j = 0; j < K; j++) {
            for (arma::uword l = 0; l < L; l++) {
                col(counter) = Theta(l, j);
                counter++;
            }
        }
    }
    if (storeW == 1 && DatObj.alphaWsFiles == 0) {
        for (arma::uword j = 0; j < K; j++) {
            arma::mat WeightsJ = Weights.slice(j);
            for (arma::uword l = 0; l < L; l++) {
                for (arma::uword o = 0; o < O; o++) {
                    for (arma::uword i = 0; i < M; i++) {
                        Index = i + M * o;
                        col(counter) = WeightsJ(l, Index);
                        counter++;
                    }
                }
            }
        }
    }

    return col;
}


arma::colvec StoreSamples(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = Para.Delta;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    arma::umat Xi = Para.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::colvec LjVec = Para.LjVec;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Xi(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = Delta(j);
        counter++;
    }    
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = LjVec(j);
        counter++;
    }
    return col;
}

arma::colvec StoreSamples(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = Para.Delta;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    arma::umat Xi = Para.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::colvec LjVec = Para.LjVec;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Xi(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = Delta(j);
        counter++;
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = LjVec(j);
        counter++;
    }
    return col;
}

arma::colvec StoreSamples(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, int numTotalPara) {

    //Set data object
    int M = DatObj.M;
    int K = DatObj.K;
    int Nu = DatObj.Nu;
    int O = DatObj.O;
    int C = DatObj.C;
    int P = DatObj.P;
    arma::Col<int> FamilyInd = DatObj.FamilyInd;

    //Set parameter objects
    arma::mat Lambda = Para.Lambda;
    arma::mat BigPhi = Para.BigPhi;
    arma::mat Sigma2 = Para.Sigma2;
    arma::mat Kappa = Para.Kappa;
    arma::colvec Delta = Para.Delta;
    arma::mat Upsilon = Para.Upsilon;
    arma::mat A = Para.A;
    arma::umat Xi = Para.Xi;
    double Rho = SpatPara.Rho;
    arma::colvec Beta = Para.Beta;
    arma::colvec LjVec = Para.LjVec;

    //Save raw samples
    int counter = 0;
    arma::uword Index;
    arma::colvec col(numTotalPara, arma::fill::zeros);
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Lambda(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword t = 0; t < Nu; t++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = BigPhi(j, t);
            counter++;
        }
    }
    if (any(FamilyInd != 3)) {
        for (arma::uword c = 0; c < (O - C); c++) {
            for (arma::uword i = 0; i < M; i++) {
                col(counter) = Sigma2(i, c);
                counter++;
            }
        }
    }
    else { counter += M * (O - C); }
    for (arma::uword i = 0; i < O; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Kappa(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j <= i; j++) {
            col(counter) = Upsilon(i, j);
            counter++;
        }
    }
    for (arma::uword i = 0; i < K; i++) {
        for (arma::uword j = 0; j < K; j++) {
            col(counter) = A(i, j);
            counter++;
        }
    }
    col(counter) = Rho;
    counter++;
    if (DatObj.P > 0) {
        for (arma::uword p = 0; p < P; p++) {
            col(counter) = Beta(p);
            counter++;
        }
    }
    else { counter += P; }
    for (arma::uword o = 0; o < O; o++) {
        for (arma::uword i = 0; i < M; i++) {
            Index = i + M * o;
            for (arma::uword j = 0; j < K; j++) {
                col(counter) = Xi(Index, j);
                counter++;
            }
        }
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = Delta(j);
        counter++;
    }
    for (arma::uword j = 0; j < K; j++) {
        col(counter) = LjVec(j);
        counter++;
    }
    return col;
}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBarInt(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgressInt);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);

  //Add percentage to submited job mode
  Rcpp::Rcout.precision(0);
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%.. ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%.. ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%.. ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%.. ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%.. ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%.. ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%.. ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%.. ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%.. ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!";

}



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBar(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  int BarLength = McmcObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0) + 1;
  Rcpp::Rcout << std::fixed << "\rBurn-in progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";
  Rcpp::Rcout << s; // for debugging
}

