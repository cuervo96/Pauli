#ifndef INTERACCION_H
#define INTERACCION_H
int Construir_Tablas(double *Tabla_VP, double *Tabla_VN, int size_tabla, double r_cut2, double s_cut2, double *deltar2, double *deltas2);
int applyPBC(double *x, double L, int N);
int PBC_delta(double *delta_X, double L);
double Calcular_E(double *Tabla_VP, double *Tabla_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, int *spin, int *particle);
double delta_E_r(double *Tabla_VP, double *Tabla_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, int target, double *x_new, int *spin, int *particle);
double delta_E_p(double *Tabla_VP, double s_cut2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, int target, double *p_new, int *spin, int *particle);
#endif
