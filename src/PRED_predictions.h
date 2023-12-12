#ifndef __predictions__
#define __predictions__

//STRUCT DEFINITIONS
struct datobjPREDtemp {
	int M;
	int K;
	int Nu;
	int O;
	int C;
	int P;
	int NNewVisits;
	arma::mat TimeDist;
	arma::uvec NewVisits;
	arma::uvec OriginalVisits;
	arma::cube Trials;
	arma::mat NewX;
	arma::Col<int> FamilyInd;
	int IT;
	int ET;
	int TempCorInd;	
	int seasonPeriod;
};
struct VAR1datobjPREDtemp {
	int M;
	int K;
	int Nu;
	int O;
	int C;
	int P;
	int NNewTime;
	arma::cube Trials;
	arma::mat NewX;
	arma::Col<int> FamilyInd;
};
struct datobjPREDspat {
	int M;
	int K;
	int Nu;
	int O;
	int C;
	int P;	
	int L;
	arma::cube Trials;
	arma::mat NewX;
	double sigma2hyparaA;
	double sigma2hyparaB;
	arma::Col<int> FamilyInd;
	int NNewLoc;
	arma::mat SpDist;
	arma::mat distOrigNew;
	arma::mat distNewNew;
	arma::umat nnIndpred;
	arma::Col<int> SeqL;
	int IS;
	int SA;
	int h;
	int SpCorInd;
	int storeW;
	int alphaWsFiles;
};
struct paraPREDtempEta {
	arma::mat Psi;
	arma::mat Upsilon;
	arma::mat Eta;
};
struct VAR1paraPREDtempEta {
	arma::mat A;
	arma::mat Upsilon;
	arma::mat Eta;
};
struct paraPREDtempY {
	arma::mat Lambda;
	arma::mat Beta;
	arma::mat Sigma2;
};
struct paraPREDspatFixedLnoCLlambda {
	arma::mat Rho;
	arma::mat Kappa;
	arma::mat Lambda;
};
struct paraPREDspatFixedLCLlambda {
	arma::mat Rho;
	arma::mat Kappa;
	arma::mat Theta;
	arma::mat Alpha;
};
struct paraPREDspatVaryLjLambda {
	arma::mat Rho;
	arma::mat Kappa;
	arma::mat LjVec;
	arma::field<arma::vec> Theta;
	arma::field<arma::mat> Alpha;
};
struct paraPREDspatY {
	arma::mat Eta;
	arma::mat Beta;
};

//PREDICTION CONVERSION FUNCTIONS
datobjPREDtemp ConvertDatObjPREDtemp(Rcpp::List DatObj_List);
VAR1datobjPREDtemp VAR1ConvertDatObjPREDtemp(Rcpp::List DatObj_List);
paraPREDtempEta ConvertParaPREDtempEta(Rcpp::List Para_List);
VAR1paraPREDtempEta VAR1ConvertParaPREDtempEta(Rcpp::List Para_List);
paraPREDtempY ConvertParaPREDtempY(Rcpp::List Para_List);
datobjPREDspat ConvertDatObjPREDspatFixedL(Rcpp::List DatObj_List);
datobjPREDspat ConvertDatObjPREDspatVaryLj(Rcpp::List DatObj_List);
paraPREDspatFixedLnoCLlambda ConvertParaPREDspatLambdaFixedLnoCL(Rcpp::List Para_List);
paraPREDspatFixedLCLlambda ConvertParaPREDspatLambdaFixedLCL(Rcpp::List Para_List);
paraPREDspatVaryLjLambda ConvertParaPREDspatLambdaVaryLj(Rcpp::List Para_List, int K, int NKeep);
paraPREDspatY ConvertParaPREDspatY(Rcpp::List Para_List);

//PREDICTION FUNCTIONS
arma::mat EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
arma::mat VAR1EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
arma::mat LambdaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
Rcpp::List AlphaKriggingFixedL(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
Rcpp::List AlphaKriggingVaryLj(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
arma::cube YKriggingTemp(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose);
arma::cube VAR1YKriggingTemp(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose);
arma::cube YKriggingSpat(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat LambdaKrig, int NKeep, bool Verbose);

//COVARIANCE FUNCTIONS
arma::mat getH(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu, int seasonPeriod);
arma::mat getInvH(double Psi, int Nu, int TempCorInd, int seasonPeriod);
arma::mat SpEXP(double Rho, arma::mat const& SpDist);
  
//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
double rbinomRcpp(double n, double p);
arma::vec rnormVecRcpp(arma::vec const& mean, arma::vec const& sd);
double pnormRcpp(double q);
double UpperpnormRcpp(double q);
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __predictions__
