int Target(N)
{
return (int) (Random() * (double)(N - 1));
}
double metropolis_r(double Energia, double *Tabla_VP, double *Tabla_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, double beta, int *aceptacion_r)
{
int i, target = Target(N);
double *x_new = (double*) malloc(3 * N * sizeof(double));
for(i = 0; i < 3 *N; i++ )
	{
	*(x_new + i) = *(x + i);
	}
double delta_r = 0.2, dE;
for (i = 0; i < 3; i++)
	{
	*(x_new + 3*target + i) = *(x_new + 3*target + i) + (Random() - 0.5) * 2 * delta_r;
	}
applyPBC(x_new, L, N);
dE = delta_E_r(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L, target, x_new);
//double E_new = Calcular_E(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x_new, p, N, rij2, pij2, L);
//double E_old = Calcular_E(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L);
//dE = E_new - E_old;
//	printf("%lf\n",dE);
if(dE < 0)
	{
	Energia += dE;
	for(i = 0; i < 3; i++)
		{
		*(x + 3 * target + i) = *(x_new + 3 * target + i);
		}
	*aceptacion_r += 1;
	}
else if(Random() < exp(-1 * beta * dE))
		{
		Energia += dE;
		for(i = 0; i < 3; i++)
			{
			*(x + 3 * target + i) = *(x_new + 3 * target + i);
			}
		*aceptacion_r += 1;
		}
free(x_new);
return Energia;
}
double metropolis_p(double Energia, double *Tabla_VP, double *Tabla_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, double beta, int *aceptacion_p)
{
int i, target = Target(N);
double *p_new = (double*) malloc(3 * N * sizeof(double));
for(i = 0; i < 3 *N; i++ )
	{
	*(p_new + i) = *(p + i);
	}
double delta_p = 1, dE;
for (i = 0; i < 3; i++)
	{
	*(p_new + 3 * target + i) = *(p_new + 3 * target + i) + (Random() - 0.5) * 2 * delta_p;
	}
dE = delta_E_p(Tabla_VP, s_cut2, deltas2, x, p, N, rij2, pij2, L, target, p_new);
//double E_new = Calcular_E(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p_new, N, rij2, pij2, L);
//double E_old = Calcular_E(Tabla_VP, Tabla_VN, r_cut2, s_cut2, deltar2, deltas2, x, p, N, rij2, pij2, L);
//dE = E_new - E_old;
if(dE < 0)
	{
	Energia += dE;
	for(i = 0; i < 3; i++)
		{
		*(p + 3 * target + i) = *(p_new + 3 * target + i);
		}
	*aceptacion_p += 1;
	}
else if(Random() < exp(-beta * dE))
		{
		Energia += dE;
		for(i = 0; i < 3; i++)
			{
			*(p + 3 * target + i) = *(p_new + 3 * target + i);
			}
		*aceptacion_p += 1;
		}
free(p_new);
return Energia;
}
