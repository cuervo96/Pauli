double Random() // returns random double between 0 and 1.
{
return (double)rand()/(double)RAND_MAX;
} 
double Gaussian(double mu, double sigma)// returns double with gaussian
{                                        //distribution.
int n = 15, i;
double z = 0;
for(i = 0; i < n; i++)
	{
	z += Random();
	}
z = sqrt(12 * (double)n)*(z/n - 0.5); // standard normal distribution.
return z * sigma + mu;
}
double delta_ij(double *x, double *delta_r, int i, int j)
{
int k;
for(k = 0; k < 3; k ++)
	{ 
	*(delta_r + k) = *(x + 3 * i + k) - *(x + 3 * j + k);
	}
return 0;
}
