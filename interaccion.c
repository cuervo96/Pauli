int Construir_Tablas(double *Tabla_VP, double *Tabla_VN, int size_tabla, double r_cut2, double s_cut2, double *deltar2, double *deltas2)
{
int i;
*deltar2 = (double)r_cut2 / (double) (size_tabla + 1); 
*deltas2 = (double)s_cut2 / (double) (size_tabla + 1);
double rij2 = *deltar2, sij2 = *deltas2;
double D = 0.0059, VP0 = 34.32, VN0= 18.263, p1 = 6.2, p2 = 3.0, r1 = 1.7456, r2 = 1.7324, alpha = 1.2, d = 3.350; 
for(i = 0; i < size_tabla; i++)
	{
	*(Tabla_VP + i) = VP0 * D * (exp(-0.5 * sij2) - exp( -0.5 * s_cut2));
	*(Tabla_VN + i) = VN0 * (pow(r1 * r1 / rij2, p1/2.0) - pow(r2 * r2 / rij2,p2/2.0)) * 1.0 / (1.0 + exp(alpha * (sqrt(rij2) - d))) - VN0 * (pow(r1 * r1 / r_cut2, p1/2.0) - pow(r2 * r2 / r_cut2,p2/2.0)) * 1.0 / (1.0 + exp(alpha * (sqrt(r_cut2) - d)));  
	rij2 += *deltar2;
	sij2 += *deltas2;
	}
return 0;
}
int applyPBC(double *x, double L, int N)
{
int i;
for(i = 0; i < 3 * N; i++)
        {
        if(*(x + i) > L)
                *(x + i) -= L;
        if(*(x + i) < 0)
                *(x + i) += L;
        }
return 0;
} 
int PBC_delta(double *delta_X, double L)
{
int i;
for(i = 0; i < 3; i++)
        {
        if(*(delta_X + i) > L/2.0)
                {
                *(delta_X + i) -= L;
                }
        else if(*(delta_X + i) <= -L/2.0)
                {
                *(delta_X + i) += L;
                }
        }
return 0;
}
double Calcular_E(double *Tabla_VP, double *Tabla_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L)
{
int i, j, k_r, k_s, k; 
double sij2, VN = 0, VP = 0, T = 0, r2 = 0, p2 = 0, q0 = 6, p0 = 2.067 * 0.0000000000000000000001;
for(i = 0; i < N - 1; i++)
	{
	for(j = i + 1; j < N; j++)	
		{r2 = 0, p2 = 0;
		delta_ij(x, rij2, i, j);
		delta_ij(p, pij2, i, j);
		//PBC_delta(rij2, L);
		for(k = 0; k < 3; k++)
			{
			r2 += *(rij2 + k) * (*(rij2 + k));
			p2 += *(pij2 + k) * (*(pij2 + k));
						
			}
//printf("r2 = %lf \n",r2);	
			if(r2 == 0)
				{
		//		printf("%lf %lf %lf %lf %lf %lf \n", *(x + 3*i), *(x + 3*i + 1), *(x + 3*i + 2),*(x + 3*j), *(x + 3*j + 1), *(x + 3*j + 2)  );
				}
		sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
	       	if(r2 < r_cut2)
			{
			k_r = (int) ((r2 - *deltar2) / *deltar2);
			if(k_r < 0)
		       		{k_r = 0;}
			//printf("%d \n", k_r );			
			VN += ((r2 - k_r * *deltar2)/ *deltar2) * (*(Tabla_VN + k_r +1)-*(Tabla_VN + k_r)) + *(Tabla_VN + k_r);		
			//printf("%lf %lf %d \n", ((r2 - k_r * *deltar2)/ *deltar2) * (*(Tabla_VN + k_r +1)-*(Tabla_VN + k_r)) + *(Tabla_VN + k_r), r2, k_r);
			}
		if(sij2 < s_cut2)
			{
			k_s = (int) ((sij2 - *deltas2) / *deltas2);
			if(k_s < 0)
				{k_s = 0;}
			VP += ((sij2 - k_s * *deltas2)/ *deltas2) * (*(Tabla_VP + k_s +1)-*(Tabla_VP + k_s)) + *(Tabla_VP + k_s);		
			}
		}
	}
for(i = 0; i < 3 *N; i++ )
	{
	T += *(p + i) * *(p + i) * 0.5;
	}
double E = T + VN + VP;	
return E;
}	
double delta_E_r(double *Tabla_VP, double *Tabla_VN, double r_cut2, double s_cut2, double *deltar2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, int target, double *x_new)
{
double dE;
int i = target, j, k_r, k_s, k; 
double sij2, VN = 0, VP = 0, r2 = 0, p2 = 0, q0 = 6, p0 = 2.067 * 0.0000000000000000000001;
for(j = 0; j < i; j++)	
	{
	delta_ij(x_new, rij2, i, j);
	delta_ij(p, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
       	if(r2 < r_cut2)
		{
		k_r = (int) ((r2 - *deltar2) / *deltar2);
		if(k_r < 0)
	       		{k_r = 0;}			
		VN += ((r2 - k_r * *deltar2) / *deltar2)* (*(Tabla_VN + k_r + 1) - *(Tabla_VN + k_r)) + *(Tabla_VN + k_r);
		}
	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP += ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for(j = i + 1; j < N; j++)	
	{
	r2 = 0;
	p2 = 0;
	delta_ij(x_new, rij2, i, j);
	delta_ij(p, pij2, i, j);
	PBC_delta(rij2, L);
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
       	if(r2 < r_cut2)
		{
		k_r = (int) ((r2 - *deltar2) / *deltar2);
		if(k_r < 0)
	       		{k_r = 0;}			
		VN += ((r2 - k_r * *deltar2) / *deltar2)* (*(Tabla_VN + k_r + 1) - *(Tabla_VN + k_r)) + *(Tabla_VN + k_r);
		}
	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP += ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for(j = 0; j < i; j++)	
	{
	delta_ij(x, rij2, i, j);
	delta_ij(p, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
       	if(r2 < r_cut2)
		{
		k_r = (int) ((r2 - *deltar2) / *deltar2);
		if(k_r < 0)
	       		{k_r = 0;}			
		VN -= ((r2 - k_r * *deltar2) / *deltar2)* (*(Tabla_VN + k_r + 1) - *(Tabla_VN + k_r)) + *(Tabla_VN + k_r);
		}
	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP -= ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for(j = i + 1; j < N; j++)	
	{
	delta_ij(x, rij2, i, j);
	delta_ij(p, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
       	if(r2 < r_cut2)
		{
		k_r = (int) ((r2 - *deltar2) / *deltar2);
		if(k_r < 0)
	       		{k_r = 0;}			
		VN -= ((r2 - k_r * *deltar2) / *deltar2)* (*(Tabla_VN + k_r + 1) - *(Tabla_VN + k_r)) + *(Tabla_VN + k_r);
		}
	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP -= ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}	
dE = VP + VN;	
return dE;
}
double delta_E_p(double *Tabla_VP, double s_cut2, double *deltas2, double *x, double *p, int N, double *rij2, double *pij2, double L, int target, double *p_new)
{
double dE;
int i = target, j, k_s, k; 
double sij2, VP = 0, T = 0, r2 = 0, p2 = 0, q0 = 6, p0 = 2.067 * 0.0000000000000000000001;
for(j = 0; j < i; j++)	
	{
	delta_ij(x, rij2, i, j);
	delta_ij(p_new, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
       	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP += ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for(j = i + 1; j < N; j++)	
	{
	delta_ij(x, rij2, i, j);
	delta_ij(p_new, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);	
	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP += ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for(j = 0; j < i; j++)	
	{
	delta_ij(x, rij2, i, j);
	delta_ij(p, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
	if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP -= ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for(j = i + 1; j < N; j++)	
	{
	delta_ij(x, rij2, i, j);
	delta_ij(p, pij2, i, j);
	PBC_delta(rij2, L);
	r2 = 0;
	p2 = 0;
	for(k = 0; k < 3; k++)
		{
		r2 += *(rij2 + k) * *(rij2 + k);
		p2 += *(pij2 + k) * *(pij2 + k);
		}
	sij2 = r2 / (q0 * q0) + p2 / (p0 * p0);
       		if(sij2 < s_cut2)
		{
		k_s = (int) ((sij2 - *deltas2) / *deltas2);
		if(k_s < 0)
			{k_s = 0;}
		VP -= ((sij2 - k_s * *deltas2) / *deltas2)* (*(Tabla_VP + k_s + 1) - *(Tabla_VP + k_s)) + *(Tabla_VP + k_s);
		}
	}
for (j = 0; j < 3 * N; j++ )
	{
	T += (*(p_new + j) * *(p_new + j) - *(p + j) * *(p + j)) * 0.5;
	}
dE = T + VP;	
return dE;
}
