#include <Rmath.h>
#include <Rcpp.h>
#include "bssa_2014_sub_header.h"

using namespace Rcpp;

//' The GMPE of BSSA 2014
//'
//' This function calculates the ground motion median values and standard deviations
//' @useDynLib RPSHA
//' @importFrom Rcpp sourceCpp
//' @param M A list of moment magnitudes for all interest of events
//' @param period The period of interest. -1 for PGV, 0 for PGA. A numeric value.
//' @param Rjb A list of Joyner-Boore distance (km); closest distance (km) to surface
//' projection of rupture plane. The length should equal to the length of moment magnitudes
//' @param Fault_Type A list of fault type indicators: 0 for unspecified fault;
//' 1 for strike-slip fault; 2 for normal fault; 3 for reverse fault.
//' The length should equal to the length of moment magnitudes
//' @param region A list of region indicators: 0 for global (incl. Taiwan); 1 for California;
//' 2 for Japan; 3 for China or Turkey; 4 for Italy. The length should equal to the length of moment magnitudes
//' @param z1 Basin depth (km): depth from the ground surface to the 1km/s shear-wave horizon.
//' -999 if unknown. A numeric value.
//' @param Vs30 Shear wave velocity averaged over top 30 m (in m/s). A numeric value.
//' @param coeffs The coefficient table of BSSA 2014. You can use the internal saved data object, bssa_2014_coeffs
//' @return A dataframe with four calculated columns, median/sigma/phi/tau. Each row is the calcualted values for each event.
//' @references Boore, D. M., Stewart, J. P., Seyhan, E., and Atkinson, G. M. (2014).
//' NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for
//' Shallow Crustal Earthquakes. Earthquake Spectra, 30(3), 1057-1085.
//' @export
// [[Rcpp::export]]
DataFrame bssa_2014_nga(NumericVector M, float period, NumericVector Rjb, NumericVector Fault_Type, NumericVector region,
                        float z1, float Vs30, NumericMatrix coeffs){

  NumericVector pre_period(107);
  pre_period[0] = -1;
  pre_period[1] = 0;
  pre_period[2] = 0.010;
  pre_period[3] = 0.020;
  pre_period[4] = 0.022;
  pre_period[5] = 0.025;
  pre_period[6] = 0.029;
  pre_period[7] = 0.030;
  pre_period[8] = 0.032;
  pre_period[9] = 0.035;
  pre_period[10] = 0.036;
  pre_period[11] = 0.040;
  pre_period[12] = 0.042;
  pre_period[13] = 0.044;
  pre_period[14] = 0.045;
  pre_period[15] = 0.046;
  pre_period[16] = 0.048;
  pre_period[17] = 0.050;
  pre_period[18] = 0.055;
  pre_period[19] = 0.060;
  pre_period[20] = 0.065;
  pre_period[21] = 0.067;
  pre_period[22] = 0.070;
  pre_period[23] = 0.075;
  pre_period[24] = 0.080;
  pre_period[25] = 0.085;
  pre_period[26] = 0.090;
  pre_period[27] = 0.095;
  pre_period[28] = 0.100;
  pre_period[29] = 0.110;
  pre_period[30] = 0.120;
  pre_period[31] = 0.130;
  pre_period[32] = 0.133;
  pre_period[33] = 0.140;
  pre_period[34] = 0.150;
  pre_period[35] = 0.160;
  pre_period[36] = 0.170;
  pre_period[37] = 0.180;
  pre_period[38] = 0.190;
  pre_period[39] = 0.200;
  pre_period[40] = 0.220;
  pre_period[41] = 0.240;
  pre_period[42] = 0.250;
  pre_period[43] = 0.260;
  pre_period[44] = 0.280;
  pre_period[45] = 0.290;
  pre_period[46] = 0.300;
  pre_period[47] = 0.320;
  pre_period[48] = 0.340;
  pre_period[49] = 0.350;
  pre_period[50] = 0.360;
  pre_period[51] = 0.380;
  pre_period[52] = 0.400;
  pre_period[53] = 0.420;
  pre_period[54] = 0.440;
  pre_period[55] = 0.450;
  pre_period[56] = 0.460;
  pre_period[57] = 0.480;
  pre_period[58] = 0.500;
  pre_period[59] = 0.550;
  pre_period[60] = 0.600;
  pre_period[61] = 0.650;
  pre_period[62] = 0.667;
  pre_period[63] = 0.700;
  pre_period[64] = 0.750;
  pre_period[65] = 0.800;
  pre_period[66] = 0.850;
  pre_period[67] = 0.900;
  pre_period[68] = 0.950;
  pre_period[69] = 1.000;
  pre_period[70] = 1.100;
  pre_period[71] = 1.200;
  pre_period[72] = 1.300;
  pre_period[73] = 1.400;
  pre_period[74] = 1.500;
  pre_period[75] = 1.600;
  pre_period[76] = 1.700;
  pre_period[77] = 1.800;
  pre_period[78] = 1.900;
  pre_period[79] = 2.000;
  pre_period[80] = 2.200;
  pre_period[81] = 2.400;
  pre_period[82] = 2.500;
  pre_period[83] = 2.600;
  pre_period[84] = 2.800;
  pre_period[85] = 3.000;
  pre_period[86] = 3.200;
  pre_period[87] = 3.400;
  pre_period[88] = 3.500;
  pre_period[89] = 3.600;
  pre_period[90] = 3.800;
  pre_period[91] = 4.000;
  pre_period[92] = 4.200;
  pre_period[93] = 4.400;
  pre_period[94] = 4.600;
  pre_period[95] = 4.800;
  pre_period[96] = 5.000;
  pre_period[97] = 5.500;
  pre_period[98] = 6.000;
  pre_period[99] = 6.500;
  pre_period[100] = 7.000;
  pre_period[101] = 7.500;
  pre_period[102] = 8.000;
  pre_period[103] = 8.500;
  pre_period[104] = 9.000;
  pre_period[105] = 9.500;
  pre_period[106] = 10.000;

  int num_events = M.length();
  NumericVector med(num_events);
  NumericVector sig(num_events);
  NumericVector phi(num_events);
  NumericVector tau(num_events);
  NumericVector result(4);
  NumericVector res_low(4);
  NumericVector res_high(4);

  int U = 0;
  int SS = 0;
  int NS = 0;
  int RS = 0;

  /* find the corresponding index to the input period */
  int ip_low = 0;
  int ip_high = pre_period.length() - 1;
  for (int j = 0; j < pre_period.length() - 1; j++) {
    if (fabs(pre_period[ip_low] - period) < 0.0001) {
      ip_high = ip_low;
      break;
    } else if (pre_period[ip_low + 1] <= period) {
      ip_low = ip_low + 1;
    }
    if (fabs(pre_period[ip_high] - period) < 0.0001) {
      ip_low = ip_high;
      break;
    } else if (pre_period[ip_high - 1] >= period) {
      ip_high = ip_high - 1;
    }
  }

  if (ip_low == ip_high) {

    /* Use pre-defined periods */
    for (int i = 0; i < num_events; i++) {
      U = 0;
      SS = 0;
      NS = 0;
      RS = 0;
      if (Fault_Type[i] == 0) {
        U = 1;
      } else if (Fault_Type[i] == 1) {
        SS = 1;
      } else if (Fault_Type[i] == 2) {
        NS = 1;
      } else if (Fault_Type[i] == 3) {
        RS = 1;
      }
      result = bssa_2014_subroutine(M[i], ip_low, Rjb[i], U, SS, NS, RS, region[i], z1, Vs30, coeffs);
      med[i] = result[0];
      sig[i] = result[1];
      phi[i] = result[2];
      tau[i] = result[3];
    }

  } else {

    /* The user defined period requires interpolation */
    for (int i = 0; i < num_events; i++) {
      U = 0;
      SS = 0;
      NS = 0;
      RS = 0;
      if (Fault_Type[i] == 0) {
        U = 1;
      } else if (Fault_Type[i] == 1) {
        SS = 1;
      } else if (Fault_Type[i] == 2) {
        NS = 1;
      } else if (Fault_Type[i] == 3) {
        RS = 1;
      }
      res_low = bssa_2014_subroutine(M[i], ip_low, Rjb[i], U, SS, NS, RS, region[i], z1, Vs30, coeffs);
      float Sa_low = res_low[0];
      float sigma_low = res_low[1];
      float phi_low = res_low[2];
      float tau_low = res_low[3];

      res_high = bssa_2014_subroutine(M[i], ip_high, Rjb[i], U, SS, NS, RS, region[i], z1, Vs30, coeffs);
      float Sa_high = res_high[0];
      float sigma_high = res_high[1];
      float phi_high = res_high[2];
      float tau_high = res_high[3];

      if (period <= pre_period[ip_low]) {
        med[i] = Sa_low;
        sig[i] = sigma_low;
        phi[i] = phi_low;
        tau[i] = tau_low;
      } else if (period >= pre_period[ip_high]) {
        med[i] = Sa_high;
        sig[i] = sigma_high;
        phi[i] = phi_high;
        tau[i] = tau_high;
      } else {
        med[i] = Sa_low + (Sa_high - Sa_low) / (pre_period[ip_high] - pre_period[ip_low]) * (period - pre_period[ip_low]);
        sig[i] = sigma_low + (sigma_high - sigma_low) / (pre_period[ip_high] - pre_period[ip_low]) * (period - pre_period[ip_low]);
        phi[i] = phi_low + (phi_high - phi_low) / (pre_period[ip_high] - pre_period[ip_low]) * (period - pre_period[ip_low]);
        tau[i] = tau_low + (tau_high - tau_low) / (pre_period[ip_high] - pre_period[ip_low]) * (period - pre_period[ip_low]);
      }
    }
  }

  /* return */
  DataFrame final_result = DataFrame::create(Named("med") = med, Named("sigma") = sig, Named("phi") = phi, Named("tau") = tau);
  return final_result;
}






