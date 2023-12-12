//
//  functions used in spatTempBFA package
//

#ifndef __spatTempBFA__
#define __spatTempBFA__

//MCMC Samplers
Rcpp::List bfaRcppFixedL(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List ParaCL_List, Rcpp::List SpatPara_List,
                       Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive);
Rcpp::List bfaRcppVaryingLjs(Rcpp::List DatObj_List, Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List SpatPara_List,
                       Rcpp::List DatAug_List, Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive);

//STRUCT DEFINITIONS
struct datobjFixedL {
    int N;
    int M;
    int Nu;
    int K;
    int L;
    int O;
    int C;
    int P;
    int GS;
    int IT;
    int IS;
    int CL;
    int ET;
    int SA;
    int seasonPeriod;
    arma::Col<int> FamilyInd;
    arma::Col<int> Indeces;
    int TempCorInd;
    int SpCorInd;
    int h;
    int alphaMethodInd;
    int spatPred;
    int storeW;
    int alphaWsFiles;
    arma::colvec YStar;
    arma::mat YStarWide;
    arma::mat SpDist;
    arma::colvec Time;
    arma::mat TimeDist;
    arma::mat EyeNu;
    arma::Col<int> SeqL;
    arma::mat X;
    arma::cube Trials;
    arma::cube Chi;  
};
struct VAR1datobjFixedL {
    int N;
    int M;
    int Nu;
    int K;
    int L;
    int O;
    int C;
    int P;
    int GS;
    int IS;
    int CL;
    int SA;
    arma::Col<int> FamilyInd;
    arma::Col<int> Indeces;
    int SpCorInd;
    int h;
    int alphaMethodInd;
    int spatPred;
    int storeW;
    int alphaWsFiles;
    arma::colvec YStar;
    arma::mat YStarWide;
    arma::mat SpDist;
    arma::Col<int> SeqL;
    arma::mat X;
    arma::cube Trials;
    arma::cube Chi;  
};
struct datobjVaryLj {
    int N;
    int M;
    int Nu;
    int K;
    int O;
    int C;
    int P;
    int GS;
    int IT;
    int IS;
    int ET;
    int SA;
    int seasonPeriod;
    int alphaSequenInd;
    arma::Col<int> FamilyInd;
    arma::Col<int> Indeces;
    int TempCorInd;
    int SpCorInd;
    int h;
    int spatPred;
    int storeW;
    arma::colvec YStar;
    arma::mat YStarWide;
    arma::mat SpDist;
    arma::colvec Time;
    arma::mat TimeDist;
    arma::mat EyeNu;
    arma::mat X;
    arma::cube Trials;
    arma::cube Chi;
};
struct VAR1datobjVaryLj {
    int N;
    int M;
    int Nu;
    int K;
    int O;
    int C;
    int P;
    int GS;
    int IS;
    int SA;
    int alphaSequenInd;
    arma::Col<int> FamilyInd;
    arma::Col<int> Indeces;
    int SpCorInd;
    int h;
    int spatPred;
    int storeW;
    arma::colvec YStar;
    arma::mat YStarWide;
    arma::mat SpDist;
    arma::mat X;
    arma::cube Trials;
    arma::cube Chi;
};
struct hypara {
    double A;
    double B;
    double SmallUpsilon;
    double A1;
    double A2;
    double APsi;
    double BPsi;
    double ARho;
    double BRho;
    double Gamma;
    double Beta;
    double Zeta;
    arma::mat Omega;
    arma::mat BigTheta;
    arma::colvec SigmaBetaInvMuBeta;
    arma::mat SigmaBetaInv;  
};
struct VAR1hypara {
    double A;
    double B;
    double SmallUpsilon;
    double A1;
    double A2;
    double ARho;
    double BRho;
    double Zeta;
    arma::mat Omega;
    arma::mat BigTheta;
    arma::colvec SigmaBetaInvMuBeta;
    arma::mat SigmaBetaInv;  
};
struct metrobj {
    double MetropPsi;
    double MetropRho;
    int AcceptanceRho;
    int AcceptancePsi;
    arma::vec OriginalTuners; 
};
struct VAR1metrobj {
    double MetropRho;
    int AcceptanceRho;
    double OriginalTuner; 
};
struct paraFixedL {
    arma::mat Sigma2;
    arma::mat Kappa;
    double Psi;
    arma::colvec Beta;
    arma::mat Upsilon;
    arma::mat UpsilonInv;
    arma::mat RootiUpsilon;
    arma::mat Lambda;
    arma::mat BigPhi;
    arma::colvec Eta;
    arma::mat HPsi;
    arma::mat RootiHPsi;
    arma::mat HPsiInv;
    arma::colvec Mean;
    arma::mat RootiKappa;
    arma::mat KappaInv;
    arma::cube Cov;
    arma::colvec XBeta;
};
struct VAR1paraFixedL {
    arma::mat Sigma2;
    arma::mat Kappa;
    arma::mat A; //the k x k matrix for VAR(1) temporal structure
    arma::colvec Beta;
    arma::mat Upsilon;
    arma::mat UpsilonInv;
    arma::mat Lambda;
    arma::mat BigPhi;
    arma::colvec Eta;
    arma::colvec Mean;
    arma::mat RootiKappa;
    arma::mat KappaInv;
    arma::cube Cov;
    arma::colvec XBeta;
};
struct paraCLfixedL {
    arma::umat Xi;
    arma::mat Theta;
    arma::colvec Delta;
    arma::colvec Tau;
    arma::cube Alpha;
    arma::cube Z;
    arma::cube Weights;
    arma::cube logWeights;
};
struct paraVaryLj {
    arma::mat Sigma2;
    arma::colvec Delta;
    arma::mat Kappa;
    double Psi;
    arma::colvec Beta;
    arma::mat Upsilon;
    arma::mat UpsilonInv;
    arma::mat RootiUpsilon;
    arma::umat Xi;
    arma::mat Lambda;
    arma::colvec Tau;
    arma::mat BigPhi;
    arma::colvec Eta;
    arma::mat HPsi;
    arma::mat RootiHPsi;
    arma::mat HPsiInv;
    arma::colvec Mean;   
    arma::mat U;
    arma::colvec LjVec;
    arma::mat RootiKappa;
    arma::mat KappaInv;
    arma::cube Cov;
    arma::colvec XBeta;
};
struct VAR1paraVaryLj {
    arma::mat Sigma2;
    arma::colvec Delta;
    arma::mat Kappa;
    arma::mat A; //the k x k matrix for VAR(1) temporal structure
    arma::colvec Beta;
    arma::mat Upsilon;
    arma::mat UpsilonInv;
    arma::umat Xi;
    arma::mat Lambda;
    arma::colvec Tau;
    arma::mat BigPhi;
    arma::colvec Eta;
    arma::colvec Mean;   
    arma::mat U;
    arma::colvec LjVec;
    arma::mat RootiKappa;
    arma::mat KappaInv;
    arma::cube Cov;
    arma::colvec XBeta;
};
struct listParaVaryLj {
    arma::field<arma::vec> Theta;
    arma::field<arma::mat> Alpha;
    arma::field<arma::mat> Weights;
};
struct spatpara1 {
    double Rho;
    arma::mat SpCovInv;
    arma::mat RootiSpCov;
    arma::colvec DwVec;
};
struct spatpara1VaryLj {
    double Rho;
    arma::mat SpCov;
    arma::mat SpCovInv;
    arma::mat RootiSpCov;
    arma::colvec DwVec;
};
struct spatpara2 {
    double Rho;
    arma::umat nnInd;
    arma::mat SpCovInv;
    double detSpCov; 
};
struct spatpara2VaryLj {
    double Rho;
    arma::mat SpCov;
    arma::umat nnInd;
    arma::mat SpCovInv;
    double detSpCov;
};
struct spatpara3 {
    double Rho;
    arma::umat nnInd;
    arma::field<arma::mat> whichJs;
    arma::mat Bs;
    arma::colvec Fs;
    arma::mat SpCovInv;
    double detSpCov;
};
struct dataug {
    int NBelow;
    int NAbove;
    arma::uvec WhichAbove;
    arma::uvec WhichBelow; 
};
struct mcmcobj {
    int NBurn;
    int NSims;
    int NThin;
    int NPilot;
    int NTotal;
    int NKeep;
    arma::vec WhichKeep;
    arma::vec WhichPilotAdapt;
    int PilotAdaptDenominator;
    arma::vec WhichBurnInProgress;
    arma::vec WhichBurnInProgressInt;
    arma::vec WhichSamplerProgress;
    arma::vec BurnInProgress;
    int BarLength;
};

//COVARIANCE FUNCTIONS
arma::mat getH(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu, int seasonPeriod);
arma::mat getRootiH(double Psi, int Nu, int TempCorInd, int seasonPeriod);
arma::mat getInvH(double Psi, int Nu, int TempCorInd, int seasonPeriod);
arma::mat SpEXP(double rho, arma::mat const& SpDist);
arma::mat GetRooti(arma::mat const& Cov);

//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::rowvec rtmvnormRcpp(arma::vec const& mean, arma::mat const& sigma, arma::vec const& lower, arma::vec const& upper);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
double pnormRcpp(double q);
double lpnormRcpp(double q);
double UpperpnormRcpp(double q);
double lUpperpnormRcpp(double q);
double qnormRcpp(double p);
double UpperqnormRcpp(double p);
double rigammaRcpp(double Alpha, double Theta);
double rgammaRcpp(double Alpha, double Theta);
arma::mat rwishRcpp(double n, arma::mat const& V);
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
double rtnormRcpp(double mean, double sd, bool Above, double x);
arma::vec pgRcpp(arma::vec const& b, arma::vec const& c);

//MCMC CONVERSION FUNCTIONS
datobjFixedL ConvertDatObjFixedL(Rcpp::List DatObj_List);
datobjVaryLj ConvertDatObjVaryLj(Rcpp::List DatObj_List);
VAR1datobjFixedL VAR1ConvertDatObjFixedL(Rcpp::List DatObj_List);
VAR1datobjVaryLj VAR1ConvertDatObjVaryLj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
metrobj ConvertMetrObj(Rcpp::List MetrObj_List);
VAR1hypara VAR1ConvertHyPara(Rcpp::List HyPara_List);
VAR1metrobj VAR1ConvertMetrObj(Rcpp::List MetrObj_List);
paraFixedL ConvertParaFixedL(Rcpp::List Para_List);
VAR1paraFixedL VAR1ConvertParaFixedL(Rcpp::List Para_List);
paraCLfixedL ConvertParaCLfixedL(Rcpp::List ParaCL_List);
paraVaryLj ConvertParaVaryLj(Rcpp::List Para_List);
VAR1paraVaryLj VAR1ConvertParaVaryLj(Rcpp::List Para_List);
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List);
dataug ConvertDatAug(Rcpp::List DatAug_List);
spatpara1 ConvertSpatPara1(Rcpp::List SpatPara_List);
spatpara1VaryLj ConvertSpatPara1VaryLj(Rcpp::List SpatPara_List);
spatpara2 ConvertSpatPara2(Rcpp::List SpatPara_List);
spatpara2VaryLj ConvertSpatPara2VaryLj(Rcpp::List SpatPara_List);
spatpara3 ConvertSpatPara3(Rcpp::List SpatPara_List);



//MCMC SAMPLER FUNCTIONS
std::pair<paraFixedL, paraCLfixedL> SampleTheta(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL);
std::pair<paraFixedL, paraCLfixedL> SampleXi(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL);
paraCLfixedL SampleZ(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL);
paraCLfixedL SampleAlpha(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int keptIter, bool storePostPara);
paraFixedL SampleLambda(datobjFixedL DatObj, paraFixedL Para, spatpara1 SpatPara);
paraCLfixedL SampleAlpha(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int keptIter, bool storePostPara);
paraFixedL SampleLambda(datobjFixedL DatObj, paraFixedL Para, spatpara2 SpatPara);
paraCLfixedL SampleAlpha(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int keptIter, bool storePostPara);
paraFixedL SampleKappa(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, hypara HyPara);
paraFixedL SampleKappa(datobjFixedL DatObj, paraFixedL Para, spatpara1 SpatPara, hypara HyPara);
paraFixedL SampleKappa(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, hypara HyPara);
paraFixedL SampleKappa(datobjFixedL DatObj, paraFixedL Para, spatpara2 SpatPara, hypara HyPara);
paraFixedL SampleKappa(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, hypara HyPara);
paraCLfixedL SampleDelta(datobjFixedL DatObj, paraFixedL Para, hypara HyPara, paraCLfixedL ParaCL);
paraFixedL SampleEta(datobjFixedL DatObj, paraFixedL Para, hypara HyPara);
paraFixedL SampleUpsilon(datobjFixedL DatObj, paraFixedL Para, hypara HyPara);
std::pair<paraFixedL, metrobj> SamplePsi(datobjFixedL DatObj, paraFixedL Para, hypara HyPara, metrobj MetrObj);
paraFixedL SampleSigma2(datobjFixedL DatObj, paraFixedL Para, hypara HyPara);
std::pair<datobjFixedL, paraFixedL> SampleY(datobjFixedL DatObj, paraFixedL Para, dataug DatAug);
std::pair<spatpara1, metrobj> SampleRho(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, hypara HyPara, metrobj MetrObj);
std::pair<spatpara2, metrobj> SampleRho(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, hypara HyPara, metrobj MetrObj);
std::pair<spatpara3, metrobj> SampleRho(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, hypara HyPara, metrobj MetrObj);
paraFixedL SampleBeta(datobjFixedL DatObj, paraFixedL Para, hypara HyPara);

std::pair<paraVaryLj, listParaVaryLj> SampleTheta(datobjVaryLj DatObj, paraVaryLj Para, listParaVaryLj ListPara);
std::pair<paraVaryLj, listParaVaryLj> SampleXi(datobjVaryLj DatObj, paraVaryLj Para, listParaVaryLj ListPara);
listParaVaryLj SampleAlpha(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, listParaVaryLj ListPara);
listParaVaryLj SampleAlpha(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, listParaVaryLj ListPara);
listParaVaryLj SampleAlpha(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, listParaVaryLj ListPara);
paraVaryLj SampleKappa(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, hypara HyPara, listParaVaryLj ListPara);
paraVaryLj SampleKappa(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, hypara HyPara, listParaVaryLj ListPara);
paraVaryLj SampleKappa(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, hypara HyPara, listParaVaryLj ListPara);
paraVaryLj SampleDelta(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara, listParaVaryLj ListPara);
paraVaryLj SampleEta(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara);
paraVaryLj SampleUpsilon(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara);
std::pair<paraVaryLj, metrobj> SamplePsi(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara, metrobj MetrObj);
paraVaryLj SampleSigma2(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara);
std::pair<datobjVaryLj, paraVaryLj> SampleY(datobjVaryLj DatObj, paraVaryLj Para, dataug DatAug);
std::pair<spatpara1VaryLj, metrobj> SampleRho(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, hypara HyPara, metrobj MetrObj, listParaVaryLj ListPara);
std::pair<spatpara2VaryLj, metrobj> SampleRho(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, hypara HyPara, metrobj MetrObj, listParaVaryLj ListPara);
std::pair<spatpara3, metrobj> SampleRho(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, hypara HyPara, metrobj MetrObj, listParaVaryLj ListPara);
paraVaryLj SampleBeta(datobjVaryLj DatObj, paraVaryLj Para, hypara HyPara);
paraVaryLj SampleU(datobjVaryLj DatObj, paraVaryLj Para, listParaVaryLj ListPara);


// for VAR(1) temporal structure
std::pair<VAR1paraFixedL, paraCLfixedL> SampleTheta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL);
std::pair<VAR1paraFixedL, paraCLfixedL> SampleXi(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL);
paraCLfixedL SampleZ(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL);
paraCLfixedL SampleAlpha(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int keptIter, bool storePostPara);
VAR1paraFixedL SampleLambda(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara1 SpatPara);
paraCLfixedL SampleAlpha(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int keptIter, bool storePostPara);
VAR1paraFixedL SampleLambda(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara2 SpatPara);
paraCLfixedL SampleAlpha(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int keptIter, bool storePostPara);
VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, VAR1hypara HyPara);
VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara1 SpatPara, VAR1hypara HyPara);
VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, VAR1hypara HyPara);
VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara2 SpatPara, VAR1hypara HyPara);
VAR1paraFixedL SampleKappa(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, VAR1hypara HyPara);
paraCLfixedL SampleDelta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara, paraCLfixedL ParaCL);
VAR1paraFixedL SampleEta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara);
VAR1paraFixedL SampleUpsilon(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara);
VAR1paraFixedL SampleA(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara);
VAR1paraFixedL SampleSigma2(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara);
std::pair<VAR1datobjFixedL, VAR1paraFixedL> SampleY(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, dataug DatAug);
std::pair<spatpara1, VAR1metrobj> SampleRho(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj);
std::pair<spatpara2, VAR1metrobj> SampleRho(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj);
std::pair<spatpara3, VAR1metrobj> SampleRho(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj);
VAR1paraFixedL SampleBeta(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, VAR1hypara HyPara);

std::pair<VAR1paraVaryLj, listParaVaryLj> SampleTheta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, listParaVaryLj ListPara);
std::pair<VAR1paraVaryLj, listParaVaryLj> SampleXi(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, listParaVaryLj ListPara);
listParaVaryLj SampleAlpha(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, listParaVaryLj ListPara);
listParaVaryLj SampleAlpha(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, listParaVaryLj ListPara);
listParaVaryLj SampleAlpha(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, listParaVaryLj ListPara);
VAR1paraVaryLj SampleKappa(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, VAR1hypara HyPara, listParaVaryLj ListPara);
VAR1paraVaryLj SampleKappa(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, VAR1hypara HyPara, listParaVaryLj ListPara);
VAR1paraVaryLj SampleKappa(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, VAR1hypara HyPara, listParaVaryLj ListPara);
VAR1paraVaryLj SampleDelta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara, listParaVaryLj ListPara);
VAR1paraVaryLj SampleEta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara);
VAR1paraVaryLj SampleUpsilon(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara);
VAR1paraVaryLj SampleA(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara);
VAR1paraVaryLj SampleSigma2(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara);
std::pair<VAR1datobjVaryLj, VAR1paraVaryLj> SampleY(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, dataug DatAug);
std::pair<spatpara1VaryLj, VAR1metrobj> SampleRho(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj, listParaVaryLj ListPara);
std::pair<spatpara2VaryLj, VAR1metrobj> SampleRho(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj, listParaVaryLj ListPara);
std::pair<spatpara3, VAR1metrobj> SampleRho(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, VAR1hypara HyPara, VAR1metrobj MetrObj, listParaVaryLj ListPara);
VAR1paraVaryLj SampleBeta(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, VAR1hypara HyPara);
VAR1paraVaryLj SampleU(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, listParaVaryLj ListPara);


  
//MCMC UTILITY FUNCTIONS
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive);
Rcpp::List OutputMetrObj(metrobj MetrObj);
Rcpp::List OutputMetrObj(VAR1metrobj MetrObj);
metrobj PilotAdaptation(metrobj MetrObj, datobjFixedL DatObj, mcmcobj McmcObj, int s);
metrobj PilotAdaptation(metrobj MetrObj, datobjVaryLj DatObj, mcmcobj McmcObj, int s);
VAR1metrobj PilotAdaptation(VAR1metrobj MetrObj, VAR1datobjFixedL DatObj, mcmcobj McmcObj, int s);
VAR1metrobj PilotAdaptation(VAR1metrobj MetrObj, VAR1datobjVaryLj DatObj, mcmcobj McmcObj, int s);
void SamplerProgress(int s, mcmcobj McmcObj);
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, spatpara1 SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, spatpara2 SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, spatpara3 SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjFixedL DatObj, paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjVaryLj DatObj, paraVaryLj Para, spatpara1VaryLj SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjVaryLj DatObj, paraVaryLj Para, spatpara2VaryLj SpatPara, int numTotalPara);
arma::colvec StoreSamples(datobjVaryLj DatObj, paraVaryLj Para, spatpara3 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara1 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara2 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, spatpara3 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara1 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara2 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjFixedL DatObj, VAR1paraFixedL Para, paraCLfixedL ParaCL, spatpara3 SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara1VaryLj SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara2VaryLj SpatPara, int numTotalPara);
arma::colvec StoreSamples(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para, spatpara3 SpatPara, int numTotalPara);
void UpdateBurnInBar(int s, mcmcobj McmcObj);
void UpdateBurnInBarInt(int s, mcmcobj McmcObj);

//UTILITY FUNCTIONS
arma::colvec UpdateLjVec(arma::mat const& U, arma::field<arma::mat> Weights, int K, int M, int O, arma::colvec LjVec);
arma::field<arma::mat> UpdateAlpha(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec);
std::pair<listParaVaryLj, paraVaryLj> InitializeListPara(datobjVaryLj DatObj, paraVaryLj Para);
std::pair<listParaVaryLj, VAR1paraVaryLj> InitializeListPara(VAR1datobjVaryLj DatObj, VAR1paraVaryLj Para);
arma::cube GetLogWeightsFixedL(arma::cube const& Alpha, int K, int M, int L, int O);
//arma::field<arma::mat> GetLogWeightsVaryLj(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec); no need anymore
arma::cube GetWeightsFixedL(arma::cube const& Alpha, int K, int M, int L, int O);
arma::field<arma::mat> GetWeightsVaryLj(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec);
arma::field<arma::mat> whichJsCalc(arma::umat nnInd);
arma::mat getCloseNNDmat(arma::mat matx);
arma::mat getCholRobust(arma::mat matx);
arma::mat CholInv(arma::mat const& Cov);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);

#endif //spatTempBFA
