#include <Rcpp.h>
using namespace Rcpp;

// A Rcpp implementation of simple leaf area and GPP models driven by meteorological data
// Xin et al. (2018) A steady-state approximation approach to simulate seasonal leaf dynamics of deciduous broadleaf forests via climate variables
// Xin et al. (2020) A simple time-stepping scheme to simulate leaf area index, phenology, and gross primary production across deciduous broadleaf forests in the eastern United States


// [[Rcpp::export]]
double iterator_rcpp(double Rg, double f_T, double f_VPD, double GPP_c, double LAI_c, double k, double epsilon_max) {
  int condition = 0;
  int it = 1;
  double fRg_s;
  double GPP_s;
  double LAI_s_ini = 3.0;
  double LAI_s = LAI_s_ini;
  while (condition == 0) {
    fRg_s = 1 - exp(-k*LAI_s);
    GPP_s = Rg*fRg_s*epsilon_max*f_T*f_VPD;
    LAI_s = fmin(GPP_s/GPP_c*LAI_c, LAI_c);
    if(it > 3){
      if(abs(LAI_s - LAI_s_ini) < 0.01|it > 100){
        condition = 1;
      }
    }
    it = it + 1;
    LAI_s_ini = LAI_s;
  }
  return LAI_s;
}

// [[Rcpp::export]]
List LAI_model_rcpp(NumericVector Tavg, NumericVector VPD, NumericVector Rg, double TMIN, double TMAX, double VPD_MIN, double VPD_MAX,
                    double epsilon_max, double k, double GPP_c, double LAI_c, double LAI_min, int W, int n_smooth){
  
  int N = Tavg.length();
  NumericVector f_T(N);
  NumericVector f_VPD(N);
  NumericVector f_Rg(N);
  NumericVector LAI_s(N);
  NumericVector LAI(N);
  NumericVector GPP(N);

  // Temperature modifier
  for(int n = 0; n < N; n++){
    f_T[n] = fmax(fmin((Tavg[n] - TMIN)/(TMAX - TMIN), 1.0), 0.0);
  }
  
  // VPD modifier
  for(int n = 0; n < N; n++){
    f_VPD[n] = fmax(fmin(1 - (VPD[n] - VPD_MIN)/(VPD_MAX - VPD_MIN), 1.0), 0.0);
  }
  
  // Loop over the days
  for(int n = 0; n < N; n++){
    LAI_s[n] = iterator_rcpp(Rg[n], f_T[n], f_VPD[n], GPP_c,LAI_c, k, epsilon_max);
  }
  
  // Apply rolling average
  LAI_s = LAI_s+LAI_min;
  LAI = LAI_s;
  
  // Use the half of window before and after each value to apply the rolling mean
  int H = (W-1)/2;

  // Apply moving average n_smooth times
  for(int S = 0; S < n_smooth; S++){
    NumericVector temp = LAI;
    for(int n = H; n < N-H+1; n++){
      LAI[n] = mean(temp[Rcpp::Range(n-H,n+H)]);
    }
  }
  
  f_Rg = 1 - exp(-k*LAI);
  GPP = Rg*f_Rg*epsilon_max*f_T*f_VPD;
  
  // Set the start of vector to NAs
  for(int n = 0; n < H; n++){
    LAI[n] = NA_REAL;
    GPP[n] = NA_REAL;
  }
  
  // Set the end of vector to NAs
  for(int n = N-H; n < N; n++){
    LAI[n] = NA_REAL;
    GPP[n] = NA_REAL;
  }
  
  // Return the outputs
  return List::create(
    _["f_T"] = f_T,
    _["f_VPD"] = f_VPD,
    _["f_Rg"] = f_Rg,
    _["LAI"] = LAI,
    _["GPP"] = GPP
  );
}