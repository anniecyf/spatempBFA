#include <RcppArmadillo.h>
#include "PRED_predictions.h"

datobjPREDtemp ConvertDatObjPREDtemp(Rcpp::List DatObj_List) {
	//Set objects from List
	int M = DatObj_List["M"];
	int K = DatObj_List["K"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int P = DatObj_List["P"];
	int NNewVisits = DatObj_List["NNewVisits"];
	arma::mat TimeDist = DatObj_List["TimeDist"];
	arma::uvec NewVisits = DatObj_List["NewVisits"];
	arma::uvec OriginalVisits = DatObj_List["OriginalVisits"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::mat NewX = DatObj_List["NewX"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int IT = DatObj_List["IT"];
	int ET = DatObj_List["ET"];	
	int TempCorInd = DatObj_List["TempCorInd"];
	int seasonPeriod = DatObj_List["seasonPeriod"];
	//Convert to C++ struct
	datobjPREDtemp DatObj;
	DatObj.M = M;
	DatObj.K = K;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.P = P;
	DatObj.NNewVisits = NNewVisits;
	DatObj.TimeDist = TimeDist;
	DatObj.NewVisits = NewVisits;
	DatObj.OriginalVisits = OriginalVisits;
	DatObj.Trials = Trials;
	DatObj.NewX = NewX;
	DatObj.FamilyInd = FamilyInd;
	DatObj.IT = IT;
	DatObj.ET = ET;
	DatObj.TempCorInd = TempCorInd;
	DatObj.seasonPeriod = seasonPeriod;	
	return DatObj;
}

VAR1datobjPREDtemp VAR1ConvertDatObjPREDtemp(Rcpp::List DatObj_List) {
	//Set objects from List
	int M = DatObj_List["M"];
	int K = DatObj_List["K"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int P = DatObj_List["P"];
	int NNewTime = DatObj_List["NNewTime"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::mat NewX = DatObj_List["NewX"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	//Convert to C++ struct
	VAR1datobjPREDtemp DatObj;
	DatObj.M = M;
	DatObj.K = K;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.P = P;
	DatObj.NNewTime = NNewTime;
	DatObj.Trials = Trials;
	DatObj.NewX = NewX;
	DatObj.FamilyInd = FamilyInd;
	return DatObj;
}

datobjPREDspat ConvertDatObjPREDspatFixedL(Rcpp::List DatObj_List) {
	//Set objects from List
	int M = DatObj_List["M"];
	int K = DatObj_List["K"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int P = DatObj_List["P"];
	int L = DatObj_List["L"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::mat NewX = DatObj_List["NewX"];
	double sigma2hyparaA = DatObj_List["sigma2hyparaA"];
	double sigma2hyparaB = DatObj_List["sigma2hyparaB"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int NNewLoc = DatObj_List["NNewLoc"];
	arma::mat SpDist = DatObj_List["SpDist"];
	arma::mat distOrigNew = DatObj_List["distOrigNew"];
	arma::mat distNewNew = DatObj_List["distNewNew"];
	arma::umat nnIndpred = DatObj_List["nnIndpred"];
	arma::Col<int> SeqL = DatObj_List["SeqL"];
	int IS = DatObj_List["IS"];
	int SA = DatObj_List["SA"];
	int h = DatObj_List["h"];
	int SpCorInd = DatObj_List["SpCorInd"];
	int storeW = DatObj_List["storeW"];
	int alphaWsFiles = DatObj_List["alphaWsFiles"];
	//Convert to C++ struct
	datobjPREDspat DatObj;
	DatObj.M = M;
	DatObj.K = K;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.P = P;
	DatObj.L = L;
	DatObj.Trials = Trials;
	DatObj.NewX = NewX;
	DatObj.sigma2hyparaA = sigma2hyparaA;
	DatObj.sigma2hyparaB = sigma2hyparaB;
	DatObj.FamilyInd = FamilyInd;
	DatObj.NNewLoc = NNewLoc;
	DatObj.SpDist = SpDist;
	DatObj.distOrigNew = distOrigNew;
	DatObj.distNewNew = distNewNew;
	DatObj.nnIndpred = nnIndpred;
	DatObj.SeqL = SeqL;
	DatObj.IS = IS;
	DatObj.SA = SA;
	DatObj.h = h;
	DatObj.SpCorInd = SpCorInd;
	DatObj.storeW = storeW;
	DatObj.alphaWsFiles = alphaWsFiles;
	return DatObj;
}

datobjPREDspat ConvertDatObjPREDspatVaryLj(Rcpp::List DatObj_List) {
	//Set objects from List
	int M = DatObj_List["M"];
	int K = DatObj_List["K"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int P = DatObj_List["P"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::mat NewX = DatObj_List["NewX"];
	double sigma2hyparaA = DatObj_List["sigma2hyparaA"];
	double sigma2hyparaB = DatObj_List["sigma2hyparaB"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int NNewLoc = DatObj_List["NNewLoc"];
	arma::mat SpDist = DatObj_List["SpDist"];
	arma::mat distOrigNew = DatObj_List["distOrigNew"];
	arma::mat distNewNew = DatObj_List["distNewNew"];
	arma::umat nnIndpred = DatObj_List["nnIndpred"];
	int IS = DatObj_List["IS"];
	int SA = DatObj_List["SA"];
	int h = DatObj_List["h"];
	int SpCorInd = DatObj_List["SpCorInd"];
	int storeW = DatObj_List["storeW"];
	//Convert to C++ struct
	datobjPREDspat DatObj{};
	DatObj.M = M;
	DatObj.K = K;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.P = P;
	DatObj.Trials = Trials;
	DatObj.NewX = NewX;
	DatObj.sigma2hyparaA = sigma2hyparaA;
	DatObj.sigma2hyparaB = sigma2hyparaB;
	DatObj.FamilyInd = FamilyInd;
	DatObj.NNewLoc = NNewLoc;
	DatObj.SpDist = SpDist;
	DatObj.distOrigNew = distOrigNew;
	DatObj.distNewNew = distNewNew;
	DatObj.nnIndpred = nnIndpred;
	DatObj.IS = IS;
	DatObj.SA = SA;
	DatObj.h = h;
	DatObj.SpCorInd = SpCorInd;
	DatObj.storeW = storeW;
	return DatObj;
}

paraPREDtempEta ConvertParaPREDtempEta(Rcpp::List Para_List) {
	//Set objects from List
	arma::mat Psi = Para_List["Psi"];
	arma::mat Upsilon = Para_List["Upsilon"];
	arma::mat Eta = Para_List["Eta"];
	//Convert to C++ struct
	paraPREDtempEta Para;
	Para.Psi = Psi;
	Para.Upsilon = Upsilon;
	Para.Eta = Eta;
	return Para;
}

VAR1paraPREDtempEta VAR1ConvertParaPREDtempEta(Rcpp::List Para_List) {
	//Set objects from List
	arma::mat A = Para_List["A"];
	arma::mat Upsilon = Para_List["Upsilon"];
	arma::mat Eta = Para_List["Eta"];
	//Convert to C++ struct
	VAR1paraPREDtempEta Para;
	Para.A = A;
	Para.Upsilon = Upsilon;
	Para.Eta = Eta;
	return Para;
}

paraPREDtempY ConvertParaPREDtempY(Rcpp::List Para_List) {
	//Set objects from List
	arma::mat Lambda = Para_List["Lambda"];
	arma::mat Beta = Para_List["Beta"];
	arma::mat Sigma2 = Para_List["Sigma2"];
	//Convert to C++ struct
	paraPREDtempY Para;
	Para.Lambda = Lambda;
	Para.Beta = Beta;
	Para.Sigma2 = Sigma2;
	return Para;
}

paraPREDspatFixedLnoCLlambda ConvertParaPREDspatLambdaFixedLnoCL(Rcpp::List Para_List) {
	//Set objects from List
	arma::mat Rho = Para_List["Rho"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::mat Lambda = Para_List["Lambda"];
	//Convert to C++ struct
	paraPREDspatFixedLnoCLlambda Para;
	Para.Rho = Rho;
	Para.Kappa = Kappa;
	Para.Lambda = Lambda;
	return Para;
}

paraPREDspatFixedLCLlambda ConvertParaPREDspatLambdaFixedLCL(Rcpp::List Para_List) {
	//Set objects from List
	arma::mat Rho = Para_List["Rho"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::mat Theta = Para_List["Theta"];
	arma::mat Alpha = Para_List["Alpha"];
	//Convert to C++ struct
	paraPREDspatFixedLCLlambda Para;
	Para.Rho = Rho;
	Para.Kappa = Kappa;
	Para.Theta = Theta;
	Para.Alpha = Alpha;
	return Para;
}

paraPREDspatVaryLjLambda ConvertParaPREDspatLambdaVaryLj(Rcpp::List Para_List, int K, int NKeep) {
	//Set objects from List
	arma::mat Rho = Para_List["Rho"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::mat LjVec = Para_List["LjVec"];
	Rcpp::List ThetaRcppList = Para_List["Theta"];
	Rcpp::List AlphaRcppList = Para_List["Alpha"];
	arma::field<arma::vec> Theta(K, NKeep);
	arma::field<arma::mat> Alpha(K, NKeep);	
	arma::uword index;
	for (arma::uword s = 0; s < NKeep; s++) {
		for (arma::uword j = 0; j < K; j++) {
			index = s * K + j;
			arma::vec ThetaElemVec = ThetaRcppList[index];
			Theta(j, s) = ThetaElemVec;
			arma::mat AlphaElemMat = AlphaRcppList[index];
			Alpha(j, s) = AlphaElemMat;
		}
	}
	//Convert to C++ struct
	paraPREDspatVaryLjLambda Para;
	Para.Rho = Rho;
	Para.Kappa = Kappa;
	Para.LjVec = LjVec;
	Para.Theta = Theta;
	Para.Alpha = Alpha;
	return Para;
}

paraPREDspatY ConvertParaPREDspatY(Rcpp::List Para_List) {
	//Set objects from List
	arma::mat Eta = Para_List["Eta"];
	arma::mat Beta = Para_List["Beta"];
	//Convert to C++ struct
	paraPREDspatY Para;
	Para.Eta = Eta;
	Para.Beta = Beta;
	return Para;
}
