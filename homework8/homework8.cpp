#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct EIGEN
{
	double *value = NULL;
	double *Q = NULL;
	double error;
}eigen;

eigen jacobi(double *A, int N, double error);//雅克比迭代，返回储存特征值，Q矩阵和误差的结构体
void printvector(double *A, int N, int n);//输出第n列的列向量
void printmatrix(double *A);//输出矩阵(debug使用)

int main()
{
	int N = 5, i;
	double error = 1e-4;
	double A[] = {3, 2, 5, 4, 6,
				  2, 1, 3,-7, 8,
				  5, 3, 2, 5,-4,
                  4,-7, 5, 1, 3,
                  6, 8,-4, 3, 8};
	eigen res = jacobi(A, N, error);
	for (i = 0; i < N; i++)
	{
		printf("r%d=%.8lf, x%d=", i + 1, res.value[i], i + 1);
		printvector(res.Q, N, i);
		printf("\n");
	}
}

eigen jacobi(double *A, int N, double error)
{
	//赋值
	double max, *res, *Q, s, t, c, d, temp;
	int i, j, p, q;
	eigen eig;
	res = (double *)malloc(N * N * sizeof(double)); 
	Q = (double *)malloc(N * N * sizeof(double));
	eig.value = (double *)malloc(N * sizeof(double));

	//寻找最大值，并对非对角元赋值
	max = 0, p = 0, q = 1;
	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			if (max < fabs(A[i * N + j]))
			{
				max = fabs(A[i * N + j]);
				p = i;
				q = j;
			}
			res[i * N + j] = A[i * N + j];
			res[j * N + i] = A[i * N + j];
			Q[i * N + j] = 0;
			Q[j * N + i] = 0;
		}
	}

	//对变换后的矩阵res和变换矩阵Q的对角元赋值
	for (i = 0; i < N; i++)
	{
		res[i * N + i] = A[i * N + i];
		Q[i * N + i] = 1;
	}

	do
	{
		//计算t，c，d，进行雅克比变换
		s =  (res[q * N + q] - res[p * N + p]) / 2 / res[p * N + q];
		c = -s + sqrt(s * s + 1);
		d = -s - sqrt(s * s + 1);
		t = (fabs(c) < fabs(d)) ? c : d;
		c = 1 / sqrt(t * t + 1);
		d = t / sqrt(t * t + 1);
		
		//对上三角部分进行操作，下三角储存上次的数据
		for (i = 0; i < p; i++)
		{
			res[i * N + p] = c * res[p * N + i] - d * res[q * N + i];
		}
		for (j = p + 1; j < N; j++)
		{
			res[p * N + j] = c * res[j * N + p] - d * res[q * N + j];
		}
		for (i = 0; i < p; i++)
		{
			res[i * N + q] = c * res[q * N + i] + d * res[p * N + i];
		}
		for (i = p + 1; i < q; i++)
		{
			res[i * N + q] = c * res[q * N + i] + d * res[i * N + p];
		}
		for (j = q + 1; j < N; j++)
		{
			res[q * N + j] = c * res[j * N + q] + d * res[j * N + p];
		}
		//计算对角元素
		res[p * N + p] -= res[q * N + p] * t;
		res[q * N + q] += res[q * N + p] * t;
		res[p * N + q] = 0;
		//将上三角部分计算结果复制到下三角部分
		for (i = 1; i < N; i++)
		{
			for (j = 0; j < i; j++)
			{
				res[i * N + j] = res[j * N + i];
			}
		}


		//计算Q = Q_1 * Q_2 ... Q_n
		//Q' = Q * Q_n
		for (i = 0; i < N; i++)
		{
			temp = Q[i * N + p];
			Q[i * N + p] = c * Q[i * N + p] - d * Q[i * N + q];
			Q[i * N + q] = c * Q[i * N + q] + d * temp;
		}

		//寻找最大值,计算A - diagA 的F范数
		max = 0, temp = 0, p = 0, q = 1;
		for (i = 0; i < N - 1; i++)
		{
			for (j = i + 1; j < N; j++)
			{
				if (max < fabs(res[i * N + j]))
				{
					max = fabs(res[i * N + j]);
					p = i;
					q = j;
				}
				temp += res[i * N + j] * res[i * N + j];
			}
		}
		temp = sqrt(2 * temp);//计算矩阵A - diag(A)的F范数

	} while (temp >= error);//判断误差是否满足条件

	eig.Q = Q;
	for (i = 0; i < N; i++)
	{
		eig.value[i] = res[i * N + i];
	}
	eig.error = max;
	return eig;
}

void printvector(double *A, int N, int n)
{
	int i;
	printf("(");
	for (i = 0; i < N - 1; i++)
	{
		printf("%.8lf,", A[i * N + n]);
	}
	printf("%.8lf", A[(N - 1) * N + n]);
	printf(")");
}

void printmatrix(double *A)
{
	int i, j;
	for (i = 0; i < 5; i++)
	{
		for (j = 0; j < 5; j++)
		{
			printf("%lf\t", A[i * 5 + j]);
		}
		printf("\n");
	}
	printf("\n");
}