#include<stdio.h>
#include<math.h>

double * gauss_jordan(double * A, double * b, int N);
int main()
{
	double A[] = { 31, -13,   0,   0,   0, -10,   0,  0,  0, 
	              -13,  35,  -9,   0, -11,   0,   0,  0,  0,
	                0,  -9,  31, -10,   0,   0,   0,  0,  0,
	                0,   0, -10,  79, -30,   0,   0,  0, -9,
	                0,   0,   0, -30,  57,  -7,   0, -5,  0,
	                0,   0,   0,   0,  -7,  47, -30,  0,  0,
	                0,   0,   0,   0,   0, -30,  41,  0,  0,
	                0,   0,   0,   0,  -5,   0,   0, 27, -2,
	                0,   0,   0,  -9,   0,   0,   0, -2, 29};
	double b[] = { -15, 27, -23, 0, -20, 12, -7, 7, 10 };
	int N = 9;
	int i;
	double *ans = gauss_jordan(A, b, N);
	for (i = 0; i < N; i++)
	{
		printf("x%d = %lf, ", i + 1, ans[i]);
	}
}

double * gauss_jordan(double * A, double * b, int N)
{
	int i, j, k, s;
	double max, temp;
	for (i = 0; i < N - 1; i++)
	{
		s = i;
		max = fabs(A[i * N + i]);
		for (k = i + 1; k < N; k++)
		{
			if (max < fabs(A[k * N + i]))
			{
				max = fabs(A[k * N + i]);
				s = k;
			}
		}

		if (s != i)
		{
			for (k = i; k < N; k++)
			{
				temp = A[i * N + k];
				A[i * N + k] = A[s * N + k];
				A[s * N + k] = temp;
				temp = b[i];
				b[i] = b[s];
				b[s] = temp;
			}
		}

		for (j = i + 1; j < N; j++)
		{
			A[j * N + i] /= A[i * N + i];
			for (k = i + 1; k < N; k++)
			{
				A[j * N + k] -= A[i * N + k] * A[j * N + i];
			}
			b[j] -= b[i] * A[j * N + i];
		}
	}

	for (i = N - 1; i >= 0; i--)
	{
		for (k = i + 1; k < N; k++)
		{
			b[i] -= b[k] * A[i * N + k];
		}
		b[i] /= A[i * N + i];
	}

	return b;
}