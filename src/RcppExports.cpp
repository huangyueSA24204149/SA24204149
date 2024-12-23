// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Gibbs
List Gibbs(int n_samples, int n, double a, double b);
RcppExport SEXP _SA24204149_Gibbs(SEXP n_samplesSEXP, SEXP nSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(Gibbs(n_samples, n, a, b));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _SA24204149_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// LossFunction
double LossFunction(NumericVector ytrue, NumericVector ypred, double alpha);
RcppExport SEXP _SA24204149_LossFunction(SEXP ytrueSEXP, SEXP ypredSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ytrue(ytrueSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ypred(ypredSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(LossFunction(ytrue, ypred, alpha));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA24204149_Gibbs", (DL_FUNC) &_SA24204149_Gibbs, 4},
    {"_SA24204149_rcpp_hello_world", (DL_FUNC) &_SA24204149_rcpp_hello_world, 0},
    {"_SA24204149_LossFunction", (DL_FUNC) &_SA24204149_LossFunction, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA24204149(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
