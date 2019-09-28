#ifndef METROPOLIS_H
#define METROPOLIS_H
int Target(int N);
double  metropolis_r(double Energia, double *Talba_VP, double *Talba_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, double beta, int *aceptacion_r, int *spin, int *particle);
double metropolis_p(double Energia, double *Talba_VP, double *Talba_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, double beta, int *aceptacion_p, int *spin, int *particle);
#endif
