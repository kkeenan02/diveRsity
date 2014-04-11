// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// glbWCcpp
List glbWCcpp(IntegerVector hsum, NumericMatrix af, NumericVector indtyp);
RcppExport SEXP diveRsity_glbWCcpp(SEXP hsumSEXP, SEXP afSEXP, SEXP indtypSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type hsum(hsumSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type indtyp(indtypSEXP );
        List __result = glbWCcpp(hsum, af, indtyp);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// myTab
NumericVector myTab(CharacterVector x);
RcppExport SEXP diveRsity_myTab(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP );
        NumericVector __result = myTab(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pwHCalc
List pwHCalc(NumericMatrix af, NumericVector sHarm, IntegerMatrix pw);
RcppExport SEXP diveRsity_pwHCalc(SEXP afSEXP, SEXP sHarmSEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sHarm(sHarmSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type pw(pwSEXP );
        List __result = pwHCalc(af, sHarm, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pwWCcpp
List pwWCcpp(List hsum1, NumericMatrix af1, NumericVector indtyp1, IntegerMatrix pw);
RcppExport SEXP diveRsity_pwWCcpp(SEXP hsum1SEXP, SEXP af1SEXP, SEXP indtyp1SEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type hsum1(hsum1SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type af1(af1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type indtyp1(indtyp1SEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type pw(pwSEXP );
        List __result = pwWCcpp(hsum1, af1, indtyp1, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// varFunc
List varFunc(NumericMatrix af, double sHarm);
RcppExport SEXP diveRsity_varFunc(SEXP afSEXP, SEXP sHarmSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< double >::type sHarm(sHarmSEXP );
        List __result = varFunc(af, sHarm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
