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

int main()
{
//----------DEFINICION VARIABLES--------------
int n = 8, N = n * n * n, pasos = 10000000, size_tabla = 500000, i;
double *x, *p, r_cut2 = 5.4 * 5.4, s_cut2 = 10, rho = 0.1, L = cbrt((double)N / rho), Temp = 2.0, beta = 1.0 / Temp;
x = (double*) malloc (3 * N * sizeof(double));
p = (double*) malloc (3 * N * sizeof(double));
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
double *E = (double*) malloc(pasos / 100 * sizeof(double));
char filename[255];
sprintf(filename, "test_metropolis.txt");
FILE *fp=fopen(filename, "w");
//---------------------------------------------
//--------CONDICIONES INICIALES----------------
Construir_Tablas(Tabla_VP, Tabla_VN, size_tabla, r_cut2, s_cut2, deltar2, deltas2);
set_box(x, N, L);
set_v(p, N, Temp);
Energia = Calcular_E(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L);
*E = Energia;
printf("%lf \n", Energia/ (double)N );
printf("%lf \n",L);
//---------------------------------------------
//------------------METROPOLIS-----------------

for(i = 1; i < pasos; i++)
	{
	Energia = metropolis_r(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_r);
	Energia = metropolis_p(Energia, Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, beta, aceptacion_p);
	if(i % 100)
		{
			*(E + i/100) = Energia;
		}
	printf("Progreso: %.2lf %% \r", (double)i / (double)pasos * 100);
	}
printf("aceptacion_r = %lf \n", (double)*aceptacion_r / (double)pasos);
printf("aceptacion_p = %lf \n", (double)*aceptacion_p / (double)pasos);
for(i = 0; i < pasos/100; i++)
	{
	fprintf(fp, "%lf \n", *(E + i)/ (double) N);
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
