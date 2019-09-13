#ifndef GENERAL_H
#define GENERAL_H
#include <math.h>
#include <time.h>
double Random();
double Gaussian(double mu, double sigma);
double delta_ij(double *x, double *delta_r, int i, int j);
#endif
