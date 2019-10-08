#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "general.h"
#include "inicializar.h"
#include "interaccion.h"
#include "metropolis.h"

#define PI 3.14159

int main(int argc,char *argv[])
{
//----------DEFINICION VARIABLES--------------
int n = 2*4, N = n * n * n, pasos = 50000,pretermalizacion = 2000000, Termalizacion = 20000, pasos_T = 100, size_tabla = 500000, i, j;
double Temp = 4.0, T_inicial = 4.0, T_final = 0.05, dT = (double)(T_inicial - T_final) / (double) pasos_T;
double *x, *p, r_cut2 = 5.4 * 5.4, s_cut2 = 10, rho,  beta = 1.0 / Temp;
sscanf(argv[1],"%lf", &rho);
double L = cbrt((double)N / rho);
x = (double*) malloc (3 * N * sizeof(double));
p = (double*) malloc (3 * N * sizeof(double));
int *spin, *particle; // ¡¡¡¡spin up = 1, spin down = -1 ---- proton = 1, neutron = -1!!!!
spin = (int*) malloc (N * sizeof(int));
particle = (int*) malloc (N * sizeof(int));
double *rij2, *pij2;
rij2 = (double*) malloc (3 * sizeof(double));
pij2 = (double*) malloc (3 * sizeof(double));
double *Tabla_VP, *Tabla_VN, *deltar2, *deltas2, Energia;
Tabla_VP = (double*) malloc(size_tabla * sizeof(double));
Tabla_VN = (double*) malloc(size_tabla * sizeof(double));
deltar2 = (double*) malloc(size_tabla * sizeof(double));
deltas2 = (double*) malloc(size_tabla * sizeof(double));
int *aceptacion_r = 0, *aceptacion_p = 0;
aceptacion_r = (int*) malloc(sizeof(int));
aceptacion_p = (int*) malloc(sizeof(int));
double *E = (double*) malloc(pasos_T * sizeof(double));
char filename[255];
sprintf(filename, "EvsT_rho=%.2lf.txt", rho);
FILE *fp=fopen(filename, "w");
//---------------------------------------------
//--------CONDICIONES INICIALES----------------
Construir_Tablas(Tabla_VP, Tabla_VN, size_tabla, r_cut2, s_cut2, deltar2, deltas2);
set_box(x, N, L);
set_v(p, N, Temp);
set_spin(spin, N);
set_particle(particle, N);
Energia = Calcular_E(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, spin, particle);
*E = Energia;
printf("%lf \n", Energia/ (double)N );
printf("%lf \n",L);
//---------------------------------------------
//------------------METROPOLIS-----------------
for(i = 0; i < pretermalizacion; i++)
	{
		Energia = metropolis_r(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_r, spin, particle);
		Energia = metropolis_p(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_p, spin, particle);
		printf("Progreso: %.2lf %% \r", (double)i / (double)(pretermalizacion + pasos_T * (Termalizacion + pasos) ) * 100);
	}
for(j = 0; j < pasos_T; j++)
	{ *aceptacion_r = 0;
		*aceptacion_p = 0;
		for(i = 0; i < Termalizacion; i++)
			{
				Energia = metropolis_r(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_r, spin, particle);
				Energia = metropolis_p(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_p, spin, particle);
				printf("Progreso: %.2lf %% \r", (double) (i + j * (Termalizacion + pasos) + pretermalizacion) / (double)(pretermalizacion + pasos_T * (Termalizacion + pasos) ) * 100);
			}
	for(i = 0; i < pasos; i++)
		{
			Energia = metropolis_r(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_r, spin, particle);
			Energia = metropolis_p(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_p, spin, particle);
			if(i % 500 == 0)
			*(E + j) += Energia * 500.0/ (double) pasos;
			printf("Progreso: %.2lf %% \r", (double) (i + (j + 1) * Termalizacion + j * pasos + pretermalizacion) / (double)(pretermalizacion + pasos_T * (Termalizacion + pasos) ) * 100);
		}
		//printf("aceptacion_r = %lf \n", (double)*aceptacion_r / (double) (pasos + Termalizacion));
		//printf("aceptacion_p = %lf \n", (double)*aceptacion_p / (double) (pasos + Termalizacion));
		Temp -= dT;
		beta = 1.0 / Temp;
	}
	for(j = 0; j < pasos_T; j++)
		{
			fprintf(fp,"%d %lf \n", j, *(E + j));
		}

//----------------------------------------------
fclose(fp);
free(x);
free(p);
free(E);
free(Tabla_VN);
free(Tabla_VP);
return 0;
}

#include "general.c"
#include "inicializar.c"
#include "interaccion.c"
#include "metropolis.c"
