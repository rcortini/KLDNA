#ifndef __KLDNA_E0_H__
#define __KLDNA_E0_H__

double k_q (double q);
double p_m (int m, double theta, double f1, double f2);
double bessel_Kn_prime (int n, double x);
double prefactor (int m, void *p);
double a_m (int m, void *p);
double b_m (int m, void *p);
double c_m (int m, void *p);

#endif
