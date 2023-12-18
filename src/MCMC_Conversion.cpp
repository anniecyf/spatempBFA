#include <RcppArmadillo.h>
#include "MCMC_bfaSpatTemp.h"

//Functions to convert Rcpp::List DatObj to custom C++ struct datobjFixedL or datobjVaryLj--------------------------------------------------
datobjFixedL ConvertDatObjFixedL(Rcpp::List DatObj_List) {

	//Set objects from List
	arma::mat YStarWide = DatObj_List["YStarWide"];
	arma::colvec YStar = DatObj_List["YStar"];
	arma::mat SpDist = DatObj_List["SpDist"];
	arma::mat TimeDist = DatObj_List["TimeDist"];
	arma::colvec Time = DatObj_List["Time"];
	int N = DatObj_List["N"];
	int K = DatObj_List["K"];
	int L = DatObj_List["L"];
	int M = DatObj_List["M"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int TempCorInd = DatObj_List["TempCorInd"];
	int SpCorInd = DatObj_List["SpCorInd"];
	int GS = DatObj_List["GS"];
	int IT = DatObj_List["IT"];
	int IS = DatObj_List["IS"];
	int CL = DatObj_List["CL"];
	int ET = DatObj_List["ET"];
	int SA = DatObj_List["SA"];
	int seasonPeriod = DatObj_List["seasonPeriod"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int h = DatObj_List["h"];
	int alphaMethodInd = DatObj_List["alphaMethodInd"];
	int spatPred = DatObj_List["spatPred"];
	int storeW = DatObj_List["storeW"];
	int alphaWsFiles = DatObj_List["alphaWsFiles"];
	arma::mat EyeNu = DatObj_List["EyeNu"];
	arma::Col<int> SeqL = DatObj_List["SeqL"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::cube Chi = DatObj_List["Chi"];
	arma::mat X = DatObj_List["X"];
	int P = DatObj_List["P"];
	arma::Col<int> Indeces = DatObj_List["Indeces"];

	//Convert to C++ struct
	datobjFixedL DatObj;
	DatObj.YStarWide = YStarWide;
	DatObj.YStar = YStar;
	DatObj.SpDist = SpDist;
	DatObj.TimeDist = TimeDist;
	DatObj.Time = Time;
	DatObj.N = N;
	DatObj.K = K;
	DatObj.L = L;
	DatObj.M = M;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.GS = GS;
	DatObj.IT = IT;
	DatObj.IS = IS;
	DatObj.CL = CL;
	DatObj.ET = ET;
	DatObj.SA = SA;
	DatObj.seasonPeriod = seasonPeriod;
	DatObj.TempCorInd = TempCorInd;
	DatObj.SpCorInd = SpCorInd;
	DatObj.FamilyInd = FamilyInd;
	DatObj.h = h;
	DatObj.alphaMethodInd = alphaMethodInd;
	DatObj.spatPred = spatPred;
	DatObj.storeW = storeW;
	DatObj.alphaWsFiles = alphaWsFiles;
	DatObj.EyeNu = EyeNu;
	DatObj.SeqL = SeqL;
	DatObj.Trials = Trials;
	DatObj.Chi = Chi;
	DatObj.P = P;
	DatObj.X = X;
	DatObj.Indeces = Indeces;
	return DatObj;

}

VAR1datobjFixedL VAR1ConvertDatObjFixedL(Rcpp::List DatObj_List) {

	//Set objects from List
	arma::mat YStarWide = DatObj_List["YStarWide"];
	arma::colvec YStar = DatObj_List["YStar"];
	arma::mat SpDist = DatObj_List["SpDist"];
	int N = DatObj_List["N"];
	int K = DatObj_List["K"];
	int L = DatObj_List["L"];
	int M = DatObj_List["M"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int SpCorInd = DatObj_List["SpCorInd"];
	int GS = DatObj_List["GS"];
	int IS = DatObj_List["IS"];
	int CL = DatObj_List["CL"];
	int SA = DatObj_List["SA"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int h = DatObj_List["h"];
	int alphaMethodInd = DatObj_List["alphaMethodInd"];
	int spatPred = DatObj_List["spatPred"];
	int storeW = DatObj_List["storeW"];
	int alphaWsFiles = DatObj_List["alphaWsFiles"];
	arma::Col<int> SeqL = DatObj_List["SeqL"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::cube Chi = DatObj_List["Chi"];
	arma::mat X = DatObj_List["X"];
	int P = DatObj_List["P"];
	arma::Col<int> Indeces = DatObj_List["Indeces"];

	//Convert to C++ struct
	VAR1datobjFixedL DatObj;
	DatObj.YStarWide = YStarWide;
	DatObj.YStar = YStar;
	DatObj.SpDist = SpDist;
	DatObj.N = N;
	DatObj.K = K;
	DatObj.L = L;
	DatObj.M = M;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.GS = GS;
	DatObj.IS = IS;
	DatObj.CL = CL;
	DatObj.SA = SA;
	DatObj.SpCorInd = SpCorInd;
	DatObj.FamilyInd = FamilyInd;
	DatObj.h = h;
	DatObj.alphaMethodInd = alphaMethodInd;
	DatObj.spatPred = spatPred;
	DatObj.storeW = storeW;
	DatObj.alphaWsFiles = alphaWsFiles;
	DatObj.SeqL = SeqL;
	DatObj.Trials = Trials;
	DatObj.Chi = Chi;
	DatObj.P = P;
	DatObj.X = X;
	DatObj.Indeces = Indeces;
	return DatObj;

}


datobjVaryLj ConvertDatObjVaryLj(Rcpp::List DatObj_List) {

	//Set objects from List
	arma::mat YStarWide = DatObj_List["YStarWide"];
	arma::colvec YStar = DatObj_List["YStar"];
	arma::mat SpDist = DatObj_List["SpDist"];
	arma::mat TimeDist = DatObj_List["TimeDist"];
	arma::colvec Time = DatObj_List["Time"];
	int N = DatObj_List["N"];
	int K = DatObj_List["K"];
	int M = DatObj_List["M"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int TempCorInd = DatObj_List["TempCorInd"];
	int SpCorInd = DatObj_List["SpCorInd"];
	int GS = DatObj_List["GS"];
	int IT = DatObj_List["IT"];
	int IS = DatObj_List["IS"];
	int ET = DatObj_List["ET"];
	int SA = DatObj_List["SA"];
	int seasonPeriod = DatObj_List["seasonPeriod"];
	int alphaSequenInd = DatObj_List["alphaSequenInd"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int h = DatObj_List["h"];
	int spatPred = DatObj_List["spatPred"];
	int storeW = DatObj_List["storeW"];
	arma::mat EyeNu = DatObj_List["EyeNu"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::cube Chi = DatObj_List["Chi"];
	arma::mat X = DatObj_List["X"];
	int P = DatObj_List["P"];
	arma::Col<int> Indeces = DatObj_List["Indeces"];

	//Convert to C++ struct
	datobjVaryLj DatObj;
	DatObj.YStarWide = YStarWide;
	DatObj.YStar = YStar;
	DatObj.SpDist = SpDist;
	DatObj.TimeDist = TimeDist;
	DatObj.Time = Time;
	DatObj.N = N;
	DatObj.K = K;
	DatObj.M = M;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.GS = GS;
	DatObj.IT = IT;
	DatObj.IS = IS;
	DatObj.ET = ET;
	DatObj.SA = SA;
	DatObj.seasonPeriod = seasonPeriod;
	DatObj.alphaSequenInd = alphaSequenInd;
	DatObj.TempCorInd = TempCorInd;
	DatObj.SpCorInd = SpCorInd;
	DatObj.FamilyInd = FamilyInd;
	DatObj.h = h;
	DatObj.spatPred = spatPred;
	DatObj.storeW = storeW;
	DatObj.EyeNu = EyeNu;
	DatObj.Trials = Trials;
	DatObj.Chi = Chi;
	DatObj.P = P;
	DatObj.X = X;
	DatObj.Indeces = Indeces;
	return DatObj;

}

VAR1datobjVaryLj VAR1ConvertDatObjVaryLj(Rcpp::List DatObj_List) {

	//Set objects from List
	arma::mat YStarWide = DatObj_List["YStarWide"];
	arma::colvec YStar = DatObj_List["YStar"];
	arma::mat SpDist = DatObj_List["SpDist"];
	int N = DatObj_List["N"];
	int K = DatObj_List["K"];
	int M = DatObj_List["M"];
	int Nu = DatObj_List["Nu"];
	int O = DatObj_List["O"];
	int C = DatObj_List["C"];
	int SpCorInd = DatObj_List["SpCorInd"];
	int GS = DatObj_List["GS"];
	int IS = DatObj_List["IS"];
	int SA = DatObj_List["SA"];
	int alphaSequenInd = DatObj_List["alphaSequenInd"];
	arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
	int h = DatObj_List["h"];
	int spatPred = DatObj_List["spatPred"];
	int storeW = DatObj_List["storeW"];
	arma::cube Trials = DatObj_List["Trials"];
	arma::cube Chi = DatObj_List["Chi"];
	arma::mat X = DatObj_List["X"];
	int P = DatObj_List["P"];
	arma::Col<int> Indeces = DatObj_List["Indeces"];

	//Convert to C++ struct
	VAR1datobjVaryLj DatObj;
	DatObj.YStarWide = YStarWide;
	DatObj.YStar = YStar;
	DatObj.SpDist = SpDist;
	DatObj.N = N;
	DatObj.K = K;
	DatObj.M = M;
	DatObj.Nu = Nu;
	DatObj.O = O;
	DatObj.C = C;
	DatObj.GS = GS;
	DatObj.IS = IS;
	DatObj.SA = SA;
	DatObj.alphaSequenInd = alphaSequenInd;
	DatObj.SpCorInd = SpCorInd;
	DatObj.FamilyInd = FamilyInd;
	DatObj.h = h;
	DatObj.spatPred = spatPred;
	DatObj.storeW = storeW;
	DatObj.Trials = Trials;
	DatObj.Chi = Chi;
	DatObj.P = P;
	DatObj.X = X;
	DatObj.Indeces = Indeces;
	return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  double A = HyPara_List["A"];
  double B = HyPara_List["B"];
  double SmallUpsilon = HyPara_List["SmallUpsilon"];
  double A1 = HyPara_List["A1"];
  double A2 = HyPara_List["A2"];
  double APsi = HyPara_List["APsi"];
  double BPsi = HyPara_List["BPsi"];
  double ARho = HyPara_List["ARho"];
  double BRho = HyPara_List["BRho"];
  double Gamma = HyPara_List["Gamma"];
  double Beta = HyPara_List["Beta"];
  double Zeta = HyPara_List["Zeta"];
  arma::mat Omega = HyPara_List["Omega"];
  arma::mat BigTheta = HyPara_List["BigTheta"];
  arma::colvec SigmaBetaInvMuBeta = HyPara_List["SigmaBetaInvMuBeta"];
  arma::mat SigmaBetaInv = HyPara_List["SigmaBetaInv"];
  
  //Convert to C++ struct
  hypara HyPara;
  HyPara.A = A;
  HyPara.B = B;
  HyPara.SmallUpsilon = SmallUpsilon;
  HyPara.A1 = A1;
  HyPara.A2 = A2;
  HyPara.APsi = APsi;
  HyPara.BPsi = BPsi;
  HyPara.ARho = ARho;
  HyPara.BRho = BRho;
  HyPara.Gamma = Gamma;
  HyPara.Beta = Beta;
  HyPara.Zeta = Zeta;
  HyPara.Omega = Omega;
  HyPara.BigTheta = BigTheta;
  HyPara.SigmaBetaInvMuBeta = SigmaBetaInvMuBeta;
  HyPara.SigmaBetaInv = SigmaBetaInv;
  return HyPara;

}

VAR1hypara VAR1ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  double A = HyPara_List["A"];
  double B = HyPara_List["B"];
  double SmallUpsilon = HyPara_List["SmallUpsilon"];
  double A1 = HyPara_List["A1"];
  double A2 = HyPara_List["A2"];
  double ARho = HyPara_List["ARho"];
  double BRho = HyPara_List["BRho"];
  double Zeta = HyPara_List["Zeta"];
  arma::mat V = HyPara_List["V"];
  arma::colvec Mvec = HyPara_List["Mvec"];
  arma::mat Omega = HyPara_List["Omega"];
  arma::mat BigTheta = HyPara_List["BigTheta"];
  arma::colvec SigmaBetaInvMuBeta = HyPara_List["SigmaBetaInvMuBeta"];
  arma::mat SigmaBetaInv = HyPara_List["SigmaBetaInv"];
  
  //Convert to C++ struct
  VAR1hypara HyPara;
  HyPara.A = A;
  HyPara.B = B;
  HyPara.SmallUpsilon = SmallUpsilon;
  HyPara.A1 = A1;
  HyPara.A2 = A2;
  HyPara.ARho = ARho;
  HyPara.BRho = BRho;
  HyPara.Zeta = Zeta;
  HyPara.V = V;
  HyPara.Mvec = Mvec;
  HyPara.Omega = Omega;
  HyPara.BigTheta = BigTheta;
  HyPara.SigmaBetaInvMuBeta = SigmaBetaInvMuBeta;
  HyPara.SigmaBetaInv = SigmaBetaInv;
  return HyPara;

}



//Function to convert Rcpp::List MetrObj to a custom C++ struct metrobj-----------------------------------------------
metrobj ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  double MetropPsi = MetrObj_List["MetropPsi"];
  int AcceptancePsi = MetrObj_List["AcceptancePsi"];
  double MetropRho = MetrObj_List["MetropRho"];
  int AcceptanceRho = MetrObj_List["AcceptanceRho"];
  arma::vec OriginalTuners = MetrObj_List["OriginalTuners"];
  
  //Convert to C++ struct
  metrobj MetrObj;
  MetrObj.MetropPsi = MetropPsi;
  MetrObj.AcceptancePsi = AcceptancePsi;
  MetrObj.MetropRho = MetropRho;
  MetrObj.AcceptanceRho = AcceptanceRho;
  MetrObj.OriginalTuners = OriginalTuners;
  return MetrObj;

}

VAR1metrobj VAR1ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  double MetropRho = MetrObj_List["MetropRho"];
  int AcceptanceRho = MetrObj_List["AcceptanceRho"];
  double OriginalTuner = MetrObj_List["OriginalTuner"];
  
  //Convert to C++ struct
  VAR1metrobj MetrObj;
  MetrObj.MetropRho = MetropRho;
  MetrObj.AcceptanceRho = AcceptanceRho;
  MetrObj.OriginalTuner = OriginalTuner;
  return MetrObj;

}


//Functions to convert Rcpp::List Para to a custom C++ struct paraFixedL or paraVaryLj-----------------------------------------------------
paraFixedL ConvertParaFixedL(Rcpp::List Para_List) {
	arma::mat Sigma2 = Para_List["Sigma2"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::colvec Beta = Para_List["Beta"];
	double Psi = Para_List["Psi"];
	arma::mat Upsilon = Para_List["Upsilon"];
	arma::mat UpsilonInv = Para_List["UpsilonInv"];
	arma::mat RootiUpsilon = Para_List["RootiUpsilon"];
	arma::mat Lambda = Para_List["Lambda"];
	arma::mat BigPhi = Para_List["BigPhi"];
	arma::colvec Eta = Para_List["Eta"];
	arma::mat HPsi = Para_List["HPsi"];
	arma::mat RootiHPsi = Para_List["RootiHPsi"];
	arma::mat HPsiInv = Para_List["HPsiInv"];
	arma::colvec Mean = Para_List["Mean"];
	arma::mat RootiKappa = Para_List["RootiKappa"];
	arma::mat KappaInv = Para_List["KappaInv"];
	arma::cube Cov = Para_List["Cov"];
	arma::colvec XBeta = Para_List["XBeta"];
	
	paraFixedL Para;
	Para.Sigma2 = Sigma2;
	Para.Kappa = Kappa;
	Para.Beta = Beta;
	Para.Psi = Psi;
	Para.Upsilon = Upsilon;
	Para.UpsilonInv = UpsilonInv;
	Para.RootiUpsilon = RootiUpsilon;
	Para.Lambda = Lambda;
	Para.BigPhi = BigPhi;
	Para.Eta = Eta;
	Para.HPsi = HPsi;
	Para.RootiHPsi = RootiHPsi;
	Para.HPsiInv = HPsiInv;
	Para.Mean = Mean;
	Para.RootiKappa = RootiKappa;
	Para.KappaInv = KappaInv;
	Para.Cov = Cov;
	Para.XBeta = XBeta;
	return Para;
}

VAR1paraFixedL VAR1ConvertParaFixedL(Rcpp::List Para_List) {
	arma::mat Sigma2 = Para_List["Sigma2"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::colvec Beta = Para_List["Beta"];
	arma::mat A = Para_List["A"];
	arma::mat Upsilon = Para_List["Upsilon"];
	arma::mat UpsilonInv = Para_List["UpsilonInv"];
	arma::mat Lambda = Para_List["Lambda"];
	arma::mat BigPhi = Para_List["BigPhi"];
	arma::colvec Eta = Para_List["Eta"];
	arma::colvec Mean = Para_List["Mean"];
	arma::mat RootiKappa = Para_List["RootiKappa"];
	arma::mat KappaInv = Para_List["KappaInv"];
	arma::cube Cov = Para_List["Cov"];
	arma::colvec XBeta = Para_List["XBeta"];
	
	VAR1paraFixedL Para;
	Para.Sigma2 = Sigma2;
	Para.Kappa = Kappa;
	Para.Beta = Beta;
	Para.A = A;
	Para.Upsilon = Upsilon;
	Para.UpsilonInv = UpsilonInv;
	Para.Lambda = Lambda;
	Para.BigPhi = BigPhi;
	Para.Eta = Eta;
	Para.Mean = Mean;
	Para.RootiKappa = RootiKappa;
	Para.KappaInv = KappaInv;
	Para.Cov = Cov;
	Para.XBeta = XBeta;
	return Para;
}


paraCLfixedL ConvertParaCLfixedL(Rcpp::List ParaCL_List) {	
	arma::umat Xi = ParaCL_List["Xi"];
	arma::mat Theta = ParaCL_List["Theta"];
	arma::colvec Delta = ParaCL_List["Delta"];
	arma::colvec Tau = ParaCL_List["Tau"];
	arma::cube Alpha = ParaCL_List["Alpha"];
	arma::cube Z = ParaCL_List["Z"];
	arma::cube Weights = ParaCL_List["Weights"];
	arma::cube logWeights = ParaCL_List["logWeights"];

	paraCLfixedL ParaCL;	
	ParaCL.Xi = Xi;
	ParaCL.Theta = Theta;
	ParaCL.Delta = Delta;
	ParaCL.Tau = Tau;	
	ParaCL.Alpha = Alpha;
	ParaCL.Z = Z;	
	ParaCL.Weights = Weights;
	ParaCL.logWeights = logWeights;
	
	return ParaCL;
}


paraVaryLj ConvertParaVaryLj(Rcpp::List Para_List) {

	//Set objects from List
	arma::mat Sigma2 = Para_List["Sigma2"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::colvec Delta = Para_List["Delta"];
	arma::colvec Beta = Para_List["Beta"];
	double Psi = Para_List["Psi"];
	arma::mat Upsilon = Para_List["Upsilon"];
	arma::mat UpsilonInv = Para_List["UpsilonInv"];
	arma::mat RootiUpsilon = Para_List["RootiUpsilon"];
	arma::umat Xi = Para_List["Xi"];
	arma::mat Lambda = Para_List["Lambda"];
	arma::colvec Tau = Para_List["Tau"];
	arma::mat BigPhi = Para_List["BigPhi"];
	arma::colvec Eta = Para_List["Eta"];
	arma::mat HPsi = Para_List["HPsi"];
	arma::mat RootiHPsi = Para_List["RootiHPsi"];
	arma::mat HPsiInv = Para_List["HPsiInv"];
	arma::colvec Mean = Para_List["Mean"];
	arma::mat U = Para_List["U"];
	arma::colvec LjVec = Para_List["LjVec"];
	arma::mat RootiKappa = Para_List["RootiKappa"];
	arma::mat KappaInv = Para_List["KappaInv"];
	arma::cube Cov = Para_List["Cov"];
	arma::colvec XBeta = Para_List["XBeta"];

	//Convert to C++ struct
	paraVaryLj Para;
	Para.Sigma2 = Sigma2;
	Para.Kappa = Kappa;
	Para.Delta = Delta;
	Para.Beta = Beta;
	Para.Psi = Psi;
	Para.Upsilon = Upsilon;
	Para.UpsilonInv = UpsilonInv;
	Para.RootiUpsilon = RootiUpsilon;
	Para.Xi = Xi;
	Para.Lambda = Lambda;
	Para.Tau = Tau;
	Para.BigPhi = BigPhi;
	Para.Eta = Eta;
	Para.HPsi = HPsi;
	Para.RootiHPsi = RootiHPsi;
	Para.HPsiInv = HPsiInv;
	Para.Mean = Mean;
	Para.U = U;
	Para.LjVec = LjVec;
	Para.RootiKappa = RootiKappa;
	Para.KappaInv = KappaInv;
	Para.Cov = Cov;
	Para.XBeta = XBeta;
	return Para;
}

VAR1paraVaryLj VAR1ConvertParaVaryLj(Rcpp::List Para_List) {

	//Set objects from List
	arma::mat Sigma2 = Para_List["Sigma2"];
	arma::mat Kappa = Para_List["Kappa"];
	arma::colvec Delta = Para_List["Delta"];
	arma::colvec Beta = Para_List["Beta"];
	arma::mat A = Para_List["A"];
	arma::mat Upsilon = Para_List["Upsilon"];
	arma::mat UpsilonInv = Para_List["UpsilonInv"];
	arma::umat Xi = Para_List["Xi"];
	arma::mat Lambda = Para_List["Lambda"];
	arma::colvec Tau = Para_List["Tau"];
	arma::mat BigPhi = Para_List["BigPhi"];
	arma::colvec Eta = Para_List["Eta"];
	arma::colvec Mean = Para_List["Mean"];
	arma::mat U = Para_List["U"];
	arma::colvec LjVec = Para_List["LjVec"];
	arma::mat RootiKappa = Para_List["RootiKappa"];
	arma::mat KappaInv = Para_List["KappaInv"];
	arma::cube Cov = Para_List["Cov"];
	arma::colvec XBeta = Para_List["XBeta"];

	//Convert to C++ struct
	VAR1paraVaryLj Para;
	Para.Sigma2 = Sigma2;
	Para.Kappa = Kappa;
	Para.Delta = Delta;
	Para.Beta = Beta;
	Para.A = A;
	Para.Upsilon = Upsilon;
	Para.UpsilonInv = UpsilonInv;
	Para.Xi = Xi;
	Para.Lambda = Lambda;
	Para.Tau = Tau;
	Para.BigPhi = BigPhi;
	Para.Eta = Eta;
	Para.Mean = Mean;
	Para.U = U;
	Para.LjVec = LjVec;
	Para.RootiKappa = RootiKappa;
	Para.KappaInv = KappaInv;
	Para.Cov = Cov;
	Para.XBeta = XBeta;
	return Para;
}



//Function to convert Rcpp::List SpatPara to a custom C++ struct spatPara1 -----------------------------------------------------
spatpara1 ConvertSpatPara1(Rcpp::List SpatPara_List) {

	//Set objects from List
	double Rho = SpatPara_List["Rho"];
	arma::mat SpCovInv = SpatPara_List["SpCovInv"];
	arma::mat RootiSpCov = SpatPara_List["RootiSpCov"];
	arma::colvec DwVec = SpatPara_List["DwVec"];

	//Convert to C++ struct
	spatpara1 SpatPara;
	SpatPara.Rho = Rho;
	SpatPara.SpCovInv = SpCovInv;
	SpatPara.RootiSpCov = RootiSpCov;
	SpatPara.DwVec = DwVec;
	return SpatPara;
}

//Function to convert Rcpp::List SpatPara to a custom C++ struct spatPara1 -----------------------------------------------------
spatpara1VaryLj ConvertSpatPara1VaryLj(Rcpp::List SpatPara_List) {

	//Set objects from List
	double Rho = SpatPara_List["Rho"];
	arma::mat SpCov = SpatPara_List["SpCov"];
	arma::mat SpCovInv = SpatPara_List["SpCovInv"];
	arma::mat RootiSpCov = SpatPara_List["RootiSpCov"];
	arma::colvec DwVec = SpatPara_List["DwVec"];

	//Convert to C++ struct
	spatpara1VaryLj SpatPara;
	SpatPara.Rho = Rho;
	SpatPara.SpCov = SpCov;
	SpatPara.SpCovInv = SpCovInv;
	SpatPara.RootiSpCov = RootiSpCov;
	SpatPara.DwVec = DwVec;
	return SpatPara;
}

// Function to convert Rcpp::List SpatPara to a custom C++ struct spatPara2---------------------------------------------------- -
spatpara2 ConvertSpatPara2(Rcpp::List SpatPara_List) {

	//Set objects from List
	double Rho = SpatPara_List["Rho"];
	arma::mat SpCovInv = SpatPara_List["SpCovInv"];
	double detSpCov = SpatPara_List["detSpCov"];
	arma::umat nnInd = SpatPara_List["nnInd"];

	//Convert to C++ struct
	spatpara2 SpatPara;
	SpatPara.Rho = Rho;
	SpatPara.SpCovInv = SpCovInv;
	SpatPara.detSpCov = detSpCov;
	SpatPara.nnInd = nnInd;
	return SpatPara;
}

// Function to convert Rcpp::List SpatPara to a custom C++ struct spatPara2VaryLj---------------------------------------------------- -
spatpara2VaryLj ConvertSpatPara2VaryLj(Rcpp::List SpatPara_List) {

	//Set objects from List
	double Rho = SpatPara_List["Rho"];
	arma::mat SpCov = SpatPara_List["SpCov"];
	arma::mat SpCovInv = SpatPara_List["SpCovInv"];
	double detSpCov = SpatPara_List["detSpCov"];
	arma::umat nnInd = SpatPara_List["nnInd"];

	//Convert to C++ struct
	spatpara2VaryLj SpatPara;
	SpatPara.Rho = Rho;
	SpatPara.SpCov = SpCov;
	SpatPara.SpCovInv = SpCovInv;
	SpatPara.detSpCov = detSpCov;
	SpatPara.nnInd = nnInd;
	return SpatPara;
}

// Function to convert Rcpp::List SpatPara to a custom C++ struct spatPara3---------------------------------------------------- -
spatpara3 ConvertSpatPara3(Rcpp::List SpatPara_List) {

	//Set objects from List
	double Rho = SpatPara_List["Rho"];
	arma::mat SpCovInv = SpatPara_List["SpCovInv"];
	double detSpCov = SpatPara_List["detSpCov"];
	arma::umat nnInd = SpatPara_List["nnInd"];
	arma::mat Bs = SpatPara_List["Bs"];
	arma::colvec Fs = SpatPara_List["Fs"];
	arma::field<arma::mat> whichJs = whichJsCalc(nnInd);

	//Convert to C++ struct
	spatpara3 SpatPara;
	SpatPara.Rho = Rho;
	SpatPara.SpCovInv = SpCovInv;
	SpatPara.detSpCov = detSpCov;
	SpatPara.nnInd = nnInd;
	SpatPara.Bs = Bs;
	SpatPara.Fs = Fs;
	SpatPara.whichJs = whichJs;
	return SpatPara;
}



//Function to convert Rcpp::List DatAug to a custom C++ struct dataug-----------------------------------------------------
dataug ConvertDatAug(Rcpp::List DatAug_List) {

  //Set objects from List
  int NBelow = DatAug_List["NBelow"];
  int NAbove = DatAug_List["NAbove"];
  arma::uvec WhichBelow = DatAug_List["WhichBelow"];
  arma::uvec WhichAbove = DatAug_List["WhichAbove"];

  //Convert to C++ struct
  dataug DatAug;
  DatAug.NBelow = NBelow;
  DatAug.NAbove = NAbove;
  DatAug.WhichBelow = WhichBelow;
  DatAug.WhichAbove = WhichAbove;
  return DatAug;
}



//Function to convert Rcpp::List McmcObj to a custom C++ struct mcmcmobj-----------------------------------------------------
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List) {

  //Set objects from List
  int NBurn = McmcObj_List["NBurn"];
  int NSims = McmcObj_List["NSims"];
  int NThin = McmcObj_List["NThin"];
  int NPilot = McmcObj_List["NPilot"];
  int NTotal = McmcObj_List["NTotal"];
  int NKeep = McmcObj_List["NKeep"];
  arma::vec WhichKeep = McmcObj_List["WhichKeep"];
  arma::vec WhichPilotAdapt = McmcObj_List["WhichPilotAdapt"];
  int PilotAdaptDenominator = McmcObj_List["PilotAdaptDenominator"];
  arma::vec WhichBurnInProgress = McmcObj_List["WhichBurnInProgress"];
  arma::vec WhichBurnInProgressInt = McmcObj_List["WhichBurnInProgressInt"];
  arma::vec WhichSamplerProgress = McmcObj_List["WhichSamplerProgress"];
  arma::vec BurnInProgress = McmcObj_List["BurnInProgress"];
  int BarLength = McmcObj_List["BarLength"];

  //Convert to C++ struct
  mcmcobj McmcObj;
  McmcObj.NBurn = NBurn;
  McmcObj.NSims = NSims;
  McmcObj.NThin = NThin;
  McmcObj.NPilot = NPilot;
  McmcObj.NTotal = NTotal;
  McmcObj.NKeep = NKeep;
  McmcObj.WhichKeep = WhichKeep;
  McmcObj.WhichPilotAdapt = WhichPilotAdapt;
  McmcObj.PilotAdaptDenominator = PilotAdaptDenominator;
  McmcObj.WhichBurnInProgress = WhichBurnInProgress;
  McmcObj.WhichBurnInProgressInt = WhichBurnInProgressInt;
  McmcObj.WhichSamplerProgress = WhichSamplerProgress;
  McmcObj.BurnInProgress = BurnInProgress;
  McmcObj.BarLength = BarLength;
  return McmcObj;
}



