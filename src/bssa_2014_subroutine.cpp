#include <Rmath.h>
#include <Rcpp.h>

using namespace Rcpp;

//' The subroutine of GMPE of BSSA 2014
//'
//' This is a subroutine of BSSA 2014
//' @useDynLib RPSHA
//' @importFrom Rcpp sourceCpp
//' @param M Moment magnitude, a numeric value
//' @param ip The index (starting from 0) of working period in the per-defined periods: (-1,
//' 0, 0.010, 0.020, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042,
//' 0.044, 0.045, 0.046, 0.048, 0.050, 0.055, 0.060, 0.065, 0.067, 0.070, 0.075,
//' 0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.133, 0.140, 0.150,
//' 0.160, 0.170, 0.180, 0.190, 0.200, 0.220, 0.240, 0.250, 0.260, 0.280, 0.290,
//' 0.300, 0.320, 0.340, 0.350, 0.360, 0.380, 0.400, 0.420, 0.440, 0.450, 0.460,
//' 0.480, 0.500, 0.550, 0.600, 0.650, 0.667, 0.700, 0.750, 0.800, 0.850, 0.900,
//' 0.950, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900,
//' 2.000, 2.200, 2.400, 2.500, 2.600, 2.800, 3.000, 3.200, 3.400, 3.500, 3.600,
//' 3.800, 4.000, 4.200, 4.400, 4.600, 4.800, 5.000, 5.500, 6.000, 6.500, 7.000,
//' 7.500, 8.000, 8.500, 9.000, 9.500, 10.000)
//' @param Rjb Joyner-Boore distance (km); closest distance (km) to surface
//' projection of rupture plane
//' @param U The indicator of fault type: 1 for unspecified fault, 0 for otherwise
//' @param SS The indicator of fault type: 1 for strike-slip fault, 0 for otherwise
//' @param NS The indicator of fault type: 1 for normal fault, 0 for otherwise
//' @param RS The indicator of fault type: 1 for reverse fault, 0 for otherwise
//' @param region Region indicator: 0 for global (incl. Taiwan); 1 for California;
//' 2 for Japan; 3 for China or Turkey; 4 for Italy
//' @param z1 Basin depth (km): depth from the ground surface to the 1km/s shear-wave horizon.
//' -999 if unknown.
//' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s)
//' @param coeffs The coefficient table of BSSA 2014. You can use the internal saved data object, bssa_2014_coeffs
//' @return A list of four calculated results, median/sigma/phi/tau
//' @references Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M. (2014).
//' NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for
//' Shallow Crustal Earthquakes. Earthquake Spectra, 30(3), 1057-1085.
//' @export
// [[Rcpp::export]]
NumericVector bssa_2014_subroutine(float M, int ip, float Rjb, int U, int SS, int NS, int RS, int region, float z1,
                                   float Vs30, NumericMatrix coeffs){

  NumericVector period = coeffs(_, 0);
  NumericVector e0 = coeffs(_, 1);
  NumericVector e1 = coeffs(_, 2);
  NumericVector e2 = coeffs(_, 3);
  NumericVector e3 = coeffs(_, 4);
  NumericVector e4 = coeffs(_, 5);
  NumericVector e5 = coeffs(_, 6);
  NumericVector e6 = coeffs(_, 7);
  NumericVector mh = coeffs(_, 8);
  NumericVector c1 = coeffs(_, 9);
  NumericVector c2 = coeffs(_, 10);
  NumericVector c3 = coeffs(_, 11);
  double mref = coeffs(0, 12);
  double rref = coeffs(0, 13);
  NumericVector h = coeffs(_, 14);
  NumericVector deltac3_gloCATW = coeffs(_, 15);
  NumericVector deltac3_CHTU = coeffs(_, 16);
  NumericVector deltac3_ITJA = coeffs(_, 17);
  NumericVector coff_c = coeffs(_, 18);
  NumericVector vc = coeffs(_, 19);
  double v_ref = coeffs(0, 20);
  double f1 = coeffs(0, 21);
  double f3 = coeffs(0, 22);
  NumericVector f4 = coeffs(_, 23);
  NumericVector f5 = coeffs(_, 24);
  NumericVector f6 = coeffs(_, 25);
  NumericVector f7 = coeffs(_, 26);
  NumericVector R1 = coeffs(_, 27);
  NumericVector R2 = coeffs(_, 28);
  NumericVector dphiR = coeffs(_, 29);
  NumericVector dphiV = coeffs(_, 30);
  double v1 = coeffs(0, 31);
  double v2 = coeffs(0, 32);
  NumericVector phi1 = coeffs(_, 33);
  NumericVector phi2 = coeffs(_, 34);
  NumericVector tau1 = coeffs(_, 35);
  NumericVector tau2 = coeffs(_, 36);

  /* define used variables */
  double F_E;
  NumericVector deltac3(coeffs.nrow());
  double r;
  double F_P;
  NumericVector pgar(4);
  double PGA_r;
  double sigma_r;
  double phi_r;
  double tau_r;
  double ln_Flin;
  double f2;
  double ln_Fnlin;
  double dz1;
  double mu_z1;
  double F_dz1;
  double F_S;
  double ln_Y;
  double med;
  double tau;
  double phi_M;
  double phi_MR;
  double phi_MRV;
  double sig;

  /* source effect */
  if(M <= mh[ip]){
    F_E = e0[ip] * U + e1[ip] * SS + e2[ip] * NS + e3[ip] * RS + e4[ip] * (M - mh[ip]) + e5[ip] * pow(M - mh[ip], 2);
  }else{
    F_E = e0[ip] * U + e1[ip] * SS + e2[ip] * NS + e3[ip] * RS + e6[ip] * (M - mh[ip]);
  }
  /* Rprintf("The source effect is: %f \n", F_E); */

  /* path effect */
  if(region == 0 | region == 1){
    deltac3 = deltac3_gloCATW;
  }else if(region == 3){
    deltac3 = deltac3_CHTU;
  }else if(region == 2 | region == 4){
    deltac3 = deltac3_ITJA;
  }
  r = pow(pow(Rjb, 2) + pow(h[ip], 2), 0.5);
  F_P = (c1[ip] + c2[ip] * (M - mref)) * log (r / rref) + (c3[ip] + deltac3[ip]) * (r - rref);
  /*  Rprintf("The path effect is: %f \n", F_P); */


  /* site effect */
  if((Vs30 != v_ref) | (ip != 1)){
    pgar = bssa_2014_subroutine(M, 1, Rjb, U, SS, NS, RS, region, z1, v_ref, coeffs);
    PGA_r = pgar[0];
    sigma_r = pgar[1];
    phi_r = pgar[2];
    tau_r = pgar[3];

    /* linear component */
    if(Vs30<= vc[ip]){
      ln_Flin = coff_c[ip] * log(Vs30 / v_ref);
    }else{
      ln_Flin = coff_c[ip] * log(vc[ip]/ v_ref);
    }

    /* nonlinear component */
    f2 = f4[ip] * (exp(f5[ip] * (R::fmin2(Vs30, 760) - 360)) - exp(f5[ip] * (760 - 360)));
    ln_Fnlin = f1 + f2*log((PGA_r + f3) / f3);

    /* basin component */
    if(z1 != -999){
      if(region == 1){  /* if in California */
        mu_z1 = exp(-7.15 / 4 * log((pow(Vs30, 4) + pow(570.94, 4)) / (pow(1360, 4) + pow(570.94, 4)))) /1000;
      }else if(region == 2){  /* if in Japan */
        mu_z1 = exp(-5.23 / 2 * log((pow(Vs30, 2) + pow(412.39, 2)) / (pow(1360, 2) + pow(412.39, 2)))) / 1000;
      }else{
        mu_z1 = exp(-7.15 / 4 * log((pow(Vs30, 4) + pow(570.94, 4)) / (pow(1360, 4) + pow(570.94, 4)))) / 1000;
      }
      dz1 = z1 - mu_z1;
    }else{
      dz1 = 0;
    }

    if(z1 != -999){
      if(period[ip] < 0.65){
        F_dz1 = 0;
      }else if((period[ip] >= 0.65) & (fabs(dz1) <= f7[ip] / f6[ip])){
        F_dz1 = f6[ip] * dz1;
      }else{
        F_dz1 = f7[ip];
      }
    }else{
      F_dz1 = 0;
    }

    /* final site effect */
    F_S = ln_Flin + ln_Fnlin + F_dz1;


    /* final GMM prediction */
    ln_Y = F_E + F_P + F_S;
    med = exp(ln_Y);

  }else{

    /* GMM at Vs_reference condition */
    F_S = 0;
    ln_Y = F_E + F_P + F_S;
    med = exp(ln_Y);
    PGA_r = med;

    ln_Flin = 0;

    f2 = f4[ip] * (exp(f5[ip] * (R::fmin2(Vs30, 760) - 360)) - exp(f5[ip] * (760 - 360)));
    ln_Fnlin = f1 + f2 * log((PGA_r + f3) / f3);

    if(z1 != -999){
      if(region == 1){  /* if in California */
        mu_z1 = exp(-7.15 / 4 * log((pow(Vs30, 4) + pow(570.94, 4)) / (pow(1360, 4) + pow(570.94, 4)))) / 1000;
      }else if(region == 2){  /* # if in Japan */
        mu_z1 = exp(-5.23 / 2 * log((pow(Vs30, 2) + pow(412.39, 2)) / (pow(1360, 2) + pow(412.39, 2)))) / 1000;
      }else{
        mu_z1 = exp(-7.15 / 4 * log((pow(Vs30, 4) + pow(570.94, 4))  / (pow(1360, 4) + pow(570.94, 4)))) / 1000;
      }
      dz1 = z1 - mu_z1;
    }else{
      dz1 = 0;
    }

    if(z1 != -999){
      if(period[ip] < 0.65){
        F_dz1 = 0;
      }else if((period[ip] >= 0.65) & (abs(dz1) <= f7[ip] / f6[ip])){
        F_dz1 = f6[ip] * dz1;
      }else{
        F_dz1 = f7[ip];
      }
    }else{
      F_dz1 = 0;
    }
  }
  /* Rprintf("The site effect is: %f \n", F_S); */

  /* aleatory variability */
  if(M <= 4.5){
    tau = tau1[ip];
    phi_M = phi1[ip];
  }else if((4.5 < M) & (M < 5.5)){
    tau = tau1[ip] + (tau2[ip] - tau1[ip]) * (M - 4.5);
    phi_M = phi1[ip] + (phi2[ip] - phi1[ip]) * (M - 4.5);
  }else {  /* M >= 5.5 */
    tau = tau2[ip];
    phi_M = phi2[ip];
  }

  if(Rjb <= R1[ip]){
    phi_MR = phi_M;
  }else if((R1[ip] < Rjb) & (Rjb <= R2[ip])){
    phi_MR = phi_M + dphiR[ip] * (log(Rjb / R1[ip]) / log(R2[ip] / R1[ip]));
  }else if(Rjb > R2[ip]){
    phi_MR = phi_M + dphiR[ip];
  }

  if(Vs30 >= v2){
    phi_MRV = phi_MR;
  }else if((v1 <= Vs30) & (Vs30 <= v2)){
    phi_MRV = phi_MR - dphiV[ip] * (log(v2 / Vs30) / log(v2 / v1));
  }else{  /* Vs30 <= v1 */
  phi_MRV = phi_MR - dphiV[ip];
  }

  sig = sqrt(pow(phi_MRV, 2) + pow(tau, 2));

  /* return */
  NumericVector results = NumericVector::create(Named("med", med), Named("sigma", sig), Named("phi", phi_MRV), Named("tau", tau));
  return results;
}

