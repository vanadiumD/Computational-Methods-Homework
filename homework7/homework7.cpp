#include<stdio.h>
#include<math.h>
#include<stdlib.h>

typedef struct RESULT
{
	double *x;
	int steps;
}result;

result gauss_seidel(double * A, double * b,double * x0, int N, double epsilon);
result SOR(double * A, double * b, double * x0, int N, double omega, double epsilon);

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
					0,   0,   0,  -9,   0,   0,   0, -2, 29 };
	double b[] = { -15, 27, -23, 0, -20, 12, -7, 7, 10 };
	double x0[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double epsilon = 1e-7;
	int N = 9;
	int i, t, best_steps;
	double best_omega;

	//guass_seidel方法
	result ans = gauss_seidel(A, b, x0, N, epsilon);
	for (i = 0; i < N; i++)
	{
		printf("x%d = %.10lf, ", i + 1, ans.x[i]);
	}
	printf("\n steps = %d\n", ans.steps);

	//利用guass_seidel结果初始化最佳步长
	best_steps = ans.steps;
	best_omega = 1;
	//SOR方法
	for (t = 1; t < 100; t++)
	{
		ans = SOR(A, b, x0, N, double(t) / 50, epsilon);
		printf("omega = %f\n", float(t) / 50);
		for (i = 0; i < N; i++)
		{
			printf("x%d = %.10lf, ", i + 1, ans.x[i]);
			if (ans.steps < best_steps)
			{
				best_steps = ans.steps;
				best_omega = double(t) / 50;
			}
		}
		printf("\n steps = %d\n", ans.steps);
	}
	printf("最佳松弛因子为 %lf, 对应步长为 %d", best_omega, best_steps);

}

//guass seidel迭代法
result gauss_seidel(double * A, double * b, double * x0, int N, double epsilon)
{
	result res;
	double temp, *x, sum, max;
	x = (double*)malloc(N * sizeof(double));
	int i, j, steps;
	//将x0数据拷贝到新地址中，防止覆盖信息
	for (i = 0; i < N; i++)
	{
		x[i] = x0[i];
	}
	steps = 0;
	do {
		max = 1;
		for (i = 0; i < N; i++)
		{
			sum = 0;
			for (j = 0; j < N; j++)
			{
				//排除对角元素
				if (j == i)
				{
					continue;
				}
				sum -= A[i * N + j] * x[j];
			}
			temp = x[i]; //temp 储存x[i]旧值
			x[i] = (sum + b[i]) / A[i * N + i];
			max = max < fabs(temp - x[i]) ? max : fabs(temp - x[i]);
		}

		steps++;
	} while (max >= epsilon);
	res.x = x;
	res.steps = steps;
	return res;
}

result SOR(double * A, double * b, double * x0, int N, double omega, double epsilon)
{
	result res;
	double temp, *x, sum, max;
	x = (double*)malloc(N * sizeof(double));
	int i, j, steps;
	//将x0数据拷贝到新地址中，防止覆盖信息
	for (i = 0; i < N; i++)
	{
		x[i] = x0[i];
	}
	steps = 0;
	do {
		max = 1;
		for (i = 0; i < N; i++)
		{
			sum = 0;
			for (j = 0; j < N; j++)
			{
				//排除对角元素
				if (j == i)
				{
					continue;
				}
				sum -= A[i * N + j] * x[j];
			}
			temp = x[i]; //temp 储存x[i]旧值
			x[i] = (1 - omega) * x[i] + omega * (sum + b[i]) / A[i * N + i];
			max = max < fabs(temp - x[i]) ? max : fabs(temp - x[i]);
		}
		steps++;
	} while (max >= epsilon);
	res.x = x;
	res.steps = steps;
	return res;
}

void copy(double *a, double *b, int N)
{

}