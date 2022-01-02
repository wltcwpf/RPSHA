#ifndef RPSHA_H
#define RPSHA_H

#include <Rcpp.h>
Rcpp::NumericVector bssa_2014_subroutine(float M, int ip, float Rjb, int U, int SS, int NS, int RS, int region, float z1, float Vs30, Rcpp::NumericMatrix coeffs);


#endif
