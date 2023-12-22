#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define N 500
#define PI acos(-1)
double * udf_lagrange(double *x, double *y, int num, double *xi, int Num);
double f(double x)
{
	return double(1) / (1 + x * x);
}
int main()
{
	int array_n[4] = { 5, 10, 20, 40 };
	int n, s, i;
	//生成均匀离散点
	double * xi = (double*)malloc((N + 1) * sizeof(double));
	double * yi = NULL;
	double error, max;
	for (i = 0; i <= N; i++)
	{
		xi[i] = -5 + 10 * double(i) / N;
	}

	printf("均匀结点差值\n");
	for (s = 0; s < 4; s++)
	{
		//生成并计算均匀差值结点：
		n = array_n[s];
		double * x = (double*)malloc((n + 1) * sizeof(double));
		double * y = (double*)malloc((n + 1) * sizeof(double));
		for (i = 0; i <= n; i++)
		{
			x[i] = -5 + double(10) * i / n;
			y[i] = f(x[i]);
		}
		//计算差值
		yi = udf_lagrange(x, y, n + 1, xi, N + 1);
		//计算误差
		max = 0;
		for (i = 0; i <= N; i++)
		{
			error = fabs(yi[i] - f(xi[i]));
			max = max > error ? max : error;
		}
		printf("n = %d, error = %lf\n", n, max);

		free(x);
		free(y);
	}


	printf("余弦结点差值\n");
	for (s = 0; s < 4; s++)
	{
		//生成并计算余弦差值结点：
		n = array_n[s];
		double * x = (double*)malloc((n + 1) * sizeof(double));
		double * y = (double*)malloc((n + 1) * sizeof(double));
		for (i = 0; i <= n; i++)
		{
			x[i] = -5 * cos(double(2 * i + 1) * PI / (2 * n + 2));
			y[i] = double(1) / (1 + x[i] * x[i]);
		}
		//计算差值
		yi = udf_lagrange(x, y, n + 1, xi, N + 1);
		//计算误差
		max = 0;
		for (i = 0; i <= N; i++)
		{
			error = fabs(yi[i] - f(xi[i]));
			max = max > error ? max : error;
		}
		printf("n = %d, error = %lf\n", n, max);

		free(x);
		free(y);
	}
}

//拉格朗日差值函数
double * udf_lagrange(double * x, double * y, int num, double *xi, int Num)
{
	double *yi, lk;
	if (Num < num)
	{
		printf("差值节点数大于差值点数\n");
		return NULL;
	}
	yi = (double *)malloc((Num) * sizeof(double));
	int i, j, k;
	//第一层循环求解差值xi处的yi值
	for (i = 0; i < Num; i++)
	{
		yi[i] = 0;
		//第二层循环计算所有差值结点
		for (j = 0; j < num; j++)
		{
			lk = 1;
			//第三层循环计算每个差值结点中的l_k(x)函数
			for (k = 0; k < num; k++)
			{
				if (k != j)
				{
					lk *= (xi[i] - x[k]) / (x[j] - x[k]);
				}
			}
			yi[i] += lk * y[j];
		}
	}
	return yi;
}