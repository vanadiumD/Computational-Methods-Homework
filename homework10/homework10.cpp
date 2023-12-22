#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define PI acos(-1)

double* FFT(double *A, int n);
double* discrete(double f(double), double g(double), double a, double b, int N); //将函数化为离散序列
double f(double x)
{
	return 0.7 * sin(2 * PI * 2 * x) + 2 * sin(2 * PI * 5 * x);
}
double g(double x)
{
	return 0;
}
int main()
{
	int i, N1 = 128, N2 = 256;
	double *a1, *a2, *g1, *g2;
	a1 = discrete(f, g, 0, 1, N1);
	a2 = discrete(f, g, 0, 1, N2);
	g1 = FFT(a1, N1);
	g2 = FFT(a2, N2);
	printf("向量g1第i个分量:\n");
	for (i = 0; i < N1; i++)
	{
		printf("x_%d = %lf, y_%d = %lf \n", i, g1[i], i, g1[i + N1]);
	}
	printf("向量g2第i个分量:\n");
	for (i = 0; i < N2; i++)
	{
		printf("x_%d = %lf, y_%d = %lf \n", i, g2[i], i, g2[i + N2]);
	}
}

double* FFT(double *A, int N)
{
	if (N < 2)
		return A;
	if (2 == N)
	{
		A[0] = 0.5 * (A[0] + A[1]);
		A[1] = 0.5 * (A[0] - A[1]);
		A[2] = 0.5 * (A[2] + A[3]);
		A[3] = 0.5 * (A[2] - A[3]);
		return A;
	}
	
	if (N > 2)
	{
		int m = N >> 1;
		double * even = (double*)malloc(N * sizeof(double));
		double * odd = (double*)malloc(N * sizeof(double));
		int i;
		double *p0, *p1, cos_omega, sin_omega;
		for (i = 0; i < N; i++)
		{
			even[i] = A[2 * i];
			odd[i] = A[2 * i + 1];
		}
		p0 = FFT(even, m);
		p1 = FFT(odd, m);
		for (i = 0; i < m; i++)
		{
			cos_omega = cos(2 * PI * i / N);
			sin_omega = sin(2 * PI * i / N);
			A[i] = 0.5 * (p0[i] + cos_omega * p1[i] + sin_omega * p1[i + m]);
			A[i + m] = 0.5 * (p0[i] - cos_omega * p1[i] - sin_omega * p1[i + m]);
			A[i + N] = 0.5 * (p0[i + m] + cos_omega * p1[i + m] - sin_omega * p1[i]);
			A[i + N  + m] = 0.5 * (p0[i + m] - cos_omega * p1[i + m] + sin_omega * p1[i]);
		}
		return A;
	}

}

double* discrete(double f(double), double g(double), double a, double b, int N)
{
	int i;
	double h = (b - a) / N;
	double * res = (double*)malloc(N * 2 * sizeof(double));
	for (i = 0; i < N; i++)
	{
		res[i] = f(a + h * i);
		res[i + N] = g(a + h * i);
	}
	return res;
}