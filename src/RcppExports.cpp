// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bssa_2014_nga
DataFrame bssa_2014_nga(NumericVector M, float period, NumericVector Rjb, NumericVector Fault_Type, NumericVector region, float z1, float Vs30, NumericMatrix coeffs);
RcppExport SEXP _RPSHA_bssa_2014_nga(SEXP MSEXP, SEXP periodSEXP, SEXP RjbSEXP, SEXP Fault_TypeSEXP, SEXP regionSEXP, SEXP z1SEXP, SEXP Vs30SEXP, SEXP coeffsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type M(MSEXP);
    Rcpp::traits::input_parameter< float >::type period(periodSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Rjb(RjbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Fault_Type(Fault_TypeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< float >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< float >::type Vs30(Vs30SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type coeffs(coeffsSEXP);
    rcpp_result_gen = Rcpp::wrap(bssa_2014_nga(M, period, Rjb, Fault_Type, region, z1, Vs30, coeffs));
    return rcpp_result_gen;
END_RCPP
}
// bssa_2014_subroutine
NumericVector bssa_2014_subroutine(float M, int ip, float Rjb, int U, int SS, int NS, int RS, int region, float z1, float Vs30, NumericMatrix coeffs);
RcppExport SEXP _RPSHA_bssa_2014_subroutine(SEXP MSEXP, SEXP ipSEXP, SEXP RjbSEXP, SEXP USEXP, SEXP SSSEXP, SEXP NSSEXP, SEXP RSSEXP, SEXP regionSEXP, SEXP z1SEXP, SEXP Vs30SEXP, SEXP coeffsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< float >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type ip(ipSEXP);
    Rcpp::traits::input_parameter< float >::type Rjb(RjbSEXP);
    Rcpp::traits::input_parameter< int >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< int >::type NS(NSSEXP);
    Rcpp::traits::input_parameter< int >::type RS(RSSEXP);
    Rcpp::traits::input_parameter< int >::type region(regionSEXP);
    Rcpp::traits::input_parameter< float >::type z1(z1SEXP);
    Rcpp::traits::input_parameter< float >::type Vs30(Vs30SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type coeffs(coeffsSEXP);
    rcpp_result_gen = Rcpp::wrap(bssa_2014_subroutine(M, ip, Rjb, U, SS, NS, RS, region, z1, Vs30, coeffs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RPSHA_bssa_2014_nga", (DL_FUNC) &_RPSHA_bssa_2014_nga, 8},
    {"_RPSHA_bssa_2014_subroutine", (DL_FUNC) &_RPSHA_bssa_2014_subroutine, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_RPSHA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
