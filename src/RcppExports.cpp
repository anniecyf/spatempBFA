// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// GetRooti
arma::mat GetRooti(arma::mat const& Cov);
RcppExport SEXP _spatempBFA_GetRooti(SEXP CovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Cov(CovSEXP);
    rcpp_result_gen = Rcpp::wrap(GetRooti(Cov));
    return rcpp_result_gen;
END_RCPP
}
// getRootiH
arma::mat getRootiH(double Psi, int Nu, int TempCorInd, int seasonPeriod);
RcppExport SEXP _spatempBFA_getRootiH(SEXP PsiSEXP, SEXP NuSEXP, SEXP TempCorIndSEXP, SEXP seasonPeriodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< int >::type Nu(NuSEXP);
    Rcpp::traits::input_parameter< int >::type TempCorInd(TempCorIndSEXP);
    Rcpp::traits::input_parameter< int >::type seasonPeriod(seasonPeriodSEXP);
    rcpp_result_gen = Rcpp::wrap(getRootiH(Psi, Nu, TempCorInd, seasonPeriod));
    return rcpp_result_gen;
END_RCPP
}
// getInvH
arma::mat getInvH(double Psi, int Nu, int TempCorInd, int seasonPeriod);
RcppExport SEXP _spatempBFA_getInvH(SEXP PsiSEXP, SEXP NuSEXP, SEXP TempCorIndSEXP, SEXP seasonPeriodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< int >::type Nu(NuSEXP);
    Rcpp::traits::input_parameter< int >::type TempCorInd(TempCorIndSEXP);
    Rcpp::traits::input_parameter< int >::type seasonPeriod(seasonPeriodSEXP);
    rcpp_result_gen = Rcpp::wrap(getInvH(Psi, Nu, TempCorInd, seasonPeriod));
    return rcpp_result_gen;
END_RCPP
}
// getH
arma::mat getH(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu, int seasonPeriod);
RcppExport SEXP _spatempBFA_getH(SEXP PsiSEXP, SEXP TempCorIndSEXP, SEXP TimeDistSEXP, SEXP NuSEXP, SEXP seasonPeriodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< int >::type TempCorInd(TempCorIndSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type TimeDist(TimeDistSEXP);
    Rcpp::traits::input_parameter< int >::type Nu(NuSEXP);
    Rcpp::traits::input_parameter< int >::type seasonPeriod(seasonPeriodSEXP);
    rcpp_result_gen = Rcpp::wrap(getH(Psi, TempCorInd, TimeDist, Nu, seasonPeriod));
    return rcpp_result_gen;
END_RCPP
}
// SpEXP
arma::mat SpEXP(double rho, arma::mat const& SpDist);
RcppExport SEXP _spatempBFA_SpEXP(SEXP rhoSEXP, SEXP SpDistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type SpDist(SpDistSEXP);
    rcpp_result_gen = Rcpp::wrap(SpEXP(rho, SpDist));
    return rcpp_result_gen;
END_RCPP
}
// GetLogLik
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_GetLogLik(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogLik(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// GetLogLikMean
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List);
RcppExport SEXP _spatempBFA_GetLogLikMean(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogLikMean(DatObj_List, Para_List));
    return rcpp_result_gen;
END_RCPP
}
// SamplePPD
arma::mat SamplePPD(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_SamplePPD(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(SamplePPD(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// VAR1bfaRcppFixedL
Rcpp::List VAR1bfaRcppFixedL(Rcpp::List DatObj_List, Rcpp::List HyPara_List, Rcpp::List MetrObj_List, Rcpp::List Para_List, Rcpp::List ParaCL_List, Rcpp::List SpatPara_List, Rcpp::List DatAug_List, Rcpp::List McmcObj_List, arma::mat RawSamples, bool Interactive);
RcppExport SEXP _spatempBFA_VAR1bfaRcppFixedL(SEXP DatObj_ListSEXP, SEXP HyPara_ListSEXP, SEXP MetrObj_ListSEXP, SEXP Para_ListSEXP, SEXP ParaCL_ListSEXP, SEXP SpatPara_ListSEXP, SEXP DatAug_ListSEXP, SEXP McmcObj_ListSEXP, SEXP RawSamplesSEXP, SEXP InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type HyPara_List(HyPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MetrObj_List(MetrObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ParaCL_List(ParaCL_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type SpatPara_List(SpatPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type DatAug_List(DatAug_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type McmcObj_List(McmcObj_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RawSamples(RawSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type Interactive(InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(VAR1bfaRcppFixedL(DatObj_List, HyPara_List, MetrObj_List, Para_List, ParaCL_List, SpatPara_List, DatAug_List, McmcObj_List, RawSamples, Interactive));
    return rcpp_result_gen;
END_RCPP
}
// VAR1bfaRcppVaryingLjs
Rcpp::List VAR1bfaRcppVaryingLjs(Rcpp::List DatObj_List, Rcpp::List HyPara_List, Rcpp::List MetrObj_List, Rcpp::List Para_List, Rcpp::List SpatPara_List, Rcpp::List DatAug_List, Rcpp::List McmcObj_List, arma::mat RawSamples, bool Interactive);
RcppExport SEXP _spatempBFA_VAR1bfaRcppVaryingLjs(SEXP DatObj_ListSEXP, SEXP HyPara_ListSEXP, SEXP MetrObj_ListSEXP, SEXP Para_ListSEXP, SEXP SpatPara_ListSEXP, SEXP DatAug_ListSEXP, SEXP McmcObj_ListSEXP, SEXP RawSamplesSEXP, SEXP InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type HyPara_List(HyPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MetrObj_List(MetrObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type SpatPara_List(SpatPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type DatAug_List(DatAug_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type McmcObj_List(McmcObj_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RawSamples(RawSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type Interactive(InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(VAR1bfaRcppVaryingLjs(DatObj_List, HyPara_List, MetrObj_List, Para_List, SpatPara_List, DatAug_List, McmcObj_List, RawSamples, Interactive));
    return rcpp_result_gen;
END_RCPP
}
// bfaRcppFixedL
Rcpp::List bfaRcppFixedL(Rcpp::List DatObj_List, Rcpp::List HyPara_List, Rcpp::List MetrObj_List, Rcpp::List Para_List, Rcpp::List ParaCL_List, Rcpp::List SpatPara_List, Rcpp::List DatAug_List, Rcpp::List McmcObj_List, arma::mat RawSamples, bool Interactive);
RcppExport SEXP _spatempBFA_bfaRcppFixedL(SEXP DatObj_ListSEXP, SEXP HyPara_ListSEXP, SEXP MetrObj_ListSEXP, SEXP Para_ListSEXP, SEXP ParaCL_ListSEXP, SEXP SpatPara_ListSEXP, SEXP DatAug_ListSEXP, SEXP McmcObj_ListSEXP, SEXP RawSamplesSEXP, SEXP InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type HyPara_List(HyPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MetrObj_List(MetrObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type ParaCL_List(ParaCL_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type SpatPara_List(SpatPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type DatAug_List(DatAug_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type McmcObj_List(McmcObj_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RawSamples(RawSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type Interactive(InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(bfaRcppFixedL(DatObj_List, HyPara_List, MetrObj_List, Para_List, ParaCL_List, SpatPara_List, DatAug_List, McmcObj_List, RawSamples, Interactive));
    return rcpp_result_gen;
END_RCPP
}
// bfaRcppVaryingLjs
Rcpp::List bfaRcppVaryingLjs(Rcpp::List DatObj_List, Rcpp::List HyPara_List, Rcpp::List MetrObj_List, Rcpp::List Para_List, Rcpp::List SpatPara_List, Rcpp::List DatAug_List, Rcpp::List McmcObj_List, arma::mat RawSamples, bool Interactive);
RcppExport SEXP _spatempBFA_bfaRcppVaryingLjs(SEXP DatObj_ListSEXP, SEXP HyPara_ListSEXP, SEXP MetrObj_ListSEXP, SEXP Para_ListSEXP, SEXP SpatPara_ListSEXP, SEXP DatAug_ListSEXP, SEXP McmcObj_ListSEXP, SEXP RawSamplesSEXP, SEXP InteractiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type HyPara_List(HyPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MetrObj_List(MetrObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type SpatPara_List(SpatPara_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type DatAug_List(DatAug_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type McmcObj_List(McmcObj_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type RawSamples(RawSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type Interactive(InteractiveSEXP);
    rcpp_result_gen = Rcpp::wrap(bfaRcppVaryingLjs(DatObj_List, HyPara_List, MetrObj_List, Para_List, SpatPara_List, DatAug_List, McmcObj_List, RawSamples, Interactive));
    return rcpp_result_gen;
END_RCPP
}
// EtaKrigging
arma::mat EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_EtaKrigging(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(EtaKrigging(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// VAR1EtaKrigging
arma::mat VAR1EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_VAR1EtaKrigging(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(VAR1EtaKrigging(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// LambdaKrigging
Rcpp::List LambdaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_LambdaKrigging(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(LambdaKrigging(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// AlphaKriggingFixedL
Rcpp::List AlphaKriggingFixedL(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_AlphaKriggingFixedL(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(AlphaKriggingFixedL(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// AlphaKriggingVaryLj
Rcpp::List AlphaKriggingVaryLj(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_AlphaKriggingVaryLj(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(AlphaKriggingVaryLj(DatObj_List, Para_List, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// YKriggingTemp
arma::cube YKriggingTemp(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_YKriggingTemp(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP EtaKrigSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type EtaKrig(EtaKrigSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(YKriggingTemp(DatObj_List, Para_List, EtaKrig, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// VAR1YKriggingTemp
arma::cube VAR1YKriggingTemp(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_VAR1YKriggingTemp(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP EtaKrigSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type EtaKrig(EtaKrigSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(VAR1YKriggingTemp(DatObj_List, Para_List, EtaKrig, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// YKriggingSpat
arma::cube YKriggingSpat(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat LambdaKrig, int NKeep, bool Verbose);
RcppExport SEXP _spatempBFA_YKriggingSpat(SEXP DatObj_ListSEXP, SEXP Para_ListSEXP, SEXP LambdaKrigSEXP, SEXP NKeepSEXP, SEXP VerboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type DatObj_List(DatObj_ListSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Para_List(Para_ListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LambdaKrig(LambdaKrigSEXP);
    Rcpp::traits::input_parameter< int >::type NKeep(NKeepSEXP);
    Rcpp::traits::input_parameter< bool >::type Verbose(VerboseSEXP);
    rcpp_result_gen = Rcpp::wrap(YKriggingSpat(DatObj_List, Para_List, LambdaKrig, NKeep, Verbose));
    return rcpp_result_gen;
END_RCPP
}
// UpdateLjVec
arma::colvec UpdateLjVec(arma::mat const& U, arma::field<arma::mat> Weights, int K, int M, int O, arma::colvec LjVec);
RcppExport SEXP _spatempBFA_UpdateLjVec(SEXP USEXP, SEXP WeightsSEXP, SEXP KSEXP, SEXP MSEXP, SEXP OSEXP, SEXP LjVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type Weights(WeightsSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type LjVec(LjVecSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateLjVec(U, Weights, K, M, O, LjVec));
    return rcpp_result_gen;
END_RCPP
}
// UpdateAlpha
arma::field<arma::mat> UpdateAlpha(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec);
RcppExport SEXP _spatempBFA_UpdateAlpha(SEXP AlphaSEXP, SEXP KSEXP, SEXP MSEXP, SEXP OSEXP, SEXP LjVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type LjVec(LjVecSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateAlpha(Alpha, K, M, O, LjVec));
    return rcpp_result_gen;
END_RCPP
}
// GetWeightsVaryLj
arma::field<arma::mat> GetWeightsVaryLj(arma::field<arma::mat> Alpha, int K, int M, int O, arma::colvec LjVec);
RcppExport SEXP _spatempBFA_GetWeightsVaryLj(SEXP AlphaSEXP, SEXP KSEXP, SEXP MSEXP, SEXP OSEXP, SEXP LjVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type LjVec(LjVecSEXP);
    rcpp_result_gen = Rcpp::wrap(GetWeightsVaryLj(Alpha, K, M, O, LjVec));
    return rcpp_result_gen;
END_RCPP
}
// GetWeightsFixedL
arma::cube GetWeightsFixedL(arma::cube const& Alpha, int K, int M, int L, int O);
RcppExport SEXP _spatempBFA_GetWeightsFixedL(SEXP AlphaSEXP, SEXP KSEXP, SEXP MSEXP, SEXP LSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(GetWeightsFixedL(Alpha, K, M, L, O));
    return rcpp_result_gen;
END_RCPP
}
// GetLogWeightsFixedL
arma::cube GetLogWeightsFixedL(arma::cube const& Alpha, int K, int M, int L, int O);
RcppExport SEXP _spatempBFA_GetLogWeightsFixedL(SEXP AlphaSEXP, SEXP KSEXP, SEXP MSEXP, SEXP LSEXP, SEXP OSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube const& >::type Alpha(AlphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type O(OSEXP);
    rcpp_result_gen = Rcpp::wrap(GetLogWeightsFixedL(Alpha, K, M, L, O));
    return rcpp_result_gen;
END_RCPP
}
// CholInv
arma::mat CholInv(arma::mat const& Cov);
RcppExport SEXP _spatempBFA_CholInv(SEXP CovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type Cov(CovSEXP);
    rcpp_result_gen = Rcpp::wrap(CholInv(Cov));
    return rcpp_result_gen;
END_RCPP
}
// whichJsCalc
arma::field<arma::mat> whichJsCalc(arma::umat nnInd);
RcppExport SEXP _spatempBFA_whichJsCalc(SEXP nnIndSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::umat >::type nnInd(nnIndSEXP);
    rcpp_result_gen = Rcpp::wrap(whichJsCalc(nnInd));
    return rcpp_result_gen;
END_RCPP
}
// getCloseNNDmat
arma::mat getCloseNNDmat(arma::mat matx);
RcppExport SEXP _spatempBFA_getCloseNNDmat(SEXP matxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type matx(matxSEXP);
    rcpp_result_gen = Rcpp::wrap(getCloseNNDmat(matx));
    return rcpp_result_gen;
END_RCPP
}
// getCholRobust
arma::mat getCholRobust(arma::mat matx);
RcppExport SEXP _spatempBFA_getCholRobust(SEXP matxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type matx(matxSEXP);
    rcpp_result_gen = Rcpp::wrap(getCholRobust(matx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatempBFA_GetRooti", (DL_FUNC) &_spatempBFA_GetRooti, 1},
    {"_spatempBFA_getRootiH", (DL_FUNC) &_spatempBFA_getRootiH, 4},
    {"_spatempBFA_getInvH", (DL_FUNC) &_spatempBFA_getInvH, 4},
    {"_spatempBFA_getH", (DL_FUNC) &_spatempBFA_getH, 5},
    {"_spatempBFA_SpEXP", (DL_FUNC) &_spatempBFA_SpEXP, 2},
    {"_spatempBFA_GetLogLik", (DL_FUNC) &_spatempBFA_GetLogLik, 4},
    {"_spatempBFA_GetLogLikMean", (DL_FUNC) &_spatempBFA_GetLogLikMean, 2},
    {"_spatempBFA_SamplePPD", (DL_FUNC) &_spatempBFA_SamplePPD, 4},
    {"_spatempBFA_VAR1bfaRcppFixedL", (DL_FUNC) &_spatempBFA_VAR1bfaRcppFixedL, 10},
    {"_spatempBFA_VAR1bfaRcppVaryingLjs", (DL_FUNC) &_spatempBFA_VAR1bfaRcppVaryingLjs, 9},
    {"_spatempBFA_bfaRcppFixedL", (DL_FUNC) &_spatempBFA_bfaRcppFixedL, 10},
    {"_spatempBFA_bfaRcppVaryingLjs", (DL_FUNC) &_spatempBFA_bfaRcppVaryingLjs, 9},
    {"_spatempBFA_EtaKrigging", (DL_FUNC) &_spatempBFA_EtaKrigging, 4},
    {"_spatempBFA_VAR1EtaKrigging", (DL_FUNC) &_spatempBFA_VAR1EtaKrigging, 4},
    {"_spatempBFA_LambdaKrigging", (DL_FUNC) &_spatempBFA_LambdaKrigging, 4},
    {"_spatempBFA_AlphaKriggingFixedL", (DL_FUNC) &_spatempBFA_AlphaKriggingFixedL, 4},
    {"_spatempBFA_AlphaKriggingVaryLj", (DL_FUNC) &_spatempBFA_AlphaKriggingVaryLj, 4},
    {"_spatempBFA_YKriggingTemp", (DL_FUNC) &_spatempBFA_YKriggingTemp, 5},
    {"_spatempBFA_VAR1YKriggingTemp", (DL_FUNC) &_spatempBFA_VAR1YKriggingTemp, 5},
    {"_spatempBFA_YKriggingSpat", (DL_FUNC) &_spatempBFA_YKriggingSpat, 5},
    {"_spatempBFA_UpdateLjVec", (DL_FUNC) &_spatempBFA_UpdateLjVec, 6},
    {"_spatempBFA_UpdateAlpha", (DL_FUNC) &_spatempBFA_UpdateAlpha, 5},
    {"_spatempBFA_GetWeightsVaryLj", (DL_FUNC) &_spatempBFA_GetWeightsVaryLj, 5},
    {"_spatempBFA_GetWeightsFixedL", (DL_FUNC) &_spatempBFA_GetWeightsFixedL, 5},
    {"_spatempBFA_GetLogWeightsFixedL", (DL_FUNC) &_spatempBFA_GetLogWeightsFixedL, 5},
    {"_spatempBFA_CholInv", (DL_FUNC) &_spatempBFA_CholInv, 1},
    {"_spatempBFA_whichJsCalc", (DL_FUNC) &_spatempBFA_whichJsCalc, 1},
    {"_spatempBFA_getCloseNNDmat", (DL_FUNC) &_spatempBFA_getCloseNNDmat, 1},
    {"_spatempBFA_getCholRobust", (DL_FUNC) &_spatempBFA_getCholRobust, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatempBFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
