#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

// all lengths in angstroms
#define DNA_RAD 10.
#define LD 7.
#define LB 7.
#define LH 1.7
#define G_DNA 0.17769
#define L 80.160
#define PHIDNA 1.2566
#ifdef R
# error "R defined"
#endif 
#define R 26.
// alpha0 = 10 degrees
#define ALPHA0 0.174532
#define SERIES_TOL 0.001
#define NMAX 10
// zero-point energy is F0 in kT
#define F0   1.75

// shame on me, a global variable
double phi;

// modified screening length
double k_q (double q) {
  return sqrt (1./(LD*LD) + q*q);
}

// angular modes of charge distribution
double p_m (int m, double theta, double f1, double f2) {
  if (m==0)
    return theta-1;
  else {
    double p = m%2 ? -1.0 : 1.;
    return f1*theta + p*f2*theta - cos (m*PHIDNA);
  }
}

// first derivative of the Bessel Kn functions
double bessel_Kn_prime (int n, double x) {
  return -1./2.*(gsl_sf_bessel_Kn (n-1,x) + gsl_sf_bessel_Kn (n+1,x));
}

double prefactor (int m, void *p) {
  double factor = m==0 ? 1. : 2.;
  double sign = m%2 ? -1.0 : 1.;
  double *par = (double *) p;
  double theta = par [1];
  double f1 = par [2];
  double f2 = par [3];
  double km = k_q (m*G_DNA);
  double ka = km*DNA_RAD;
  double Kmp = bessel_Kn_prime (m, ka);
  double pm = p_m (m, theta, f1, f2);
  return LB*L/(LH*LH) * factor * sign * cos (m*phi) * pm*pm / (ka * Kmp*Kmp);
}

// zero-order
double a_m (int m, void *p) {
  double k = prefactor (m,p);
  double km = k_q (m*G_DNA);
  double ka = km*DNA_RAD;
  double kR = km*R;
  return 2.*k*gsl_sf_bessel_Kn (0, kR)/ka;
}

// first-order
double b_m (int m, void *p) {
  double k = prefactor (m,p);
  double km = k_q (m*G_DNA);
  double ka = km*DNA_RAD;
  double kR = km*R;
  return 2.*k*G_DNA*DNA_RAD* m*m * gsl_sf_bessel_Kn (1, kR)/(ka*ka);
}

// second-order
double c_m (int m, void *p) {
  double k = prefactor (m,p);
  double km = k_q (m*G_DNA);
  double kR = km*R;
  return -k*L*L/(6.*R*DNA_RAD)*gsl_sf_bessel_Kn (1, kR);
}
