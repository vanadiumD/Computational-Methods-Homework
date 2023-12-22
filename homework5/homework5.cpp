#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define N 5

//定义储存计算结果的结构体
typedef struct RESULT
{
	int steps;
	double *xk;
	double *fk;
	double *error;
	double x0;
} result;

double f(double x);
double f1(double x);
result newton_method(double x0, double epsilon);
result secant_method(double x0, double x1, double epsilon);

int main()
{
	//初始化变量
	double x0[4] = { 0.1, 0.2, 0.9, 9.0 };
	double x1[4] = { -0.1, -0.2, -0.2, -0.9 };
	double x2[4] = { 0.1, 0.2, 0.9, 9.0 };
	double epsilon = 1e-8;
	int i, j;
	int m = (int)(sizeof(x0) / sizeof(double));
	int n = (int)(sizeof(x1) / sizeof(double));
	result *resn = (result*)malloc(m * sizeof(result));
	result *resc = (result*)malloc(n * sizeof(result));

	//牛顿法计算并输出结果
	printf("牛顿法\n");
	for (i = 0; i < m; i++)
	{
		resn[i] = newton_method(x0[i], epsilon);
		printf("初值 = %.14lf, 根 = %.14lf, 迭代步数 = %d\n", x0[i], resn[i].x0, resn[i].steps);
	}
	printf("\n");

	//弦截法计算并输出结果
	printf("弦截法\n");
	 n = (int)(sizeof(x1) / sizeof(double));
	for (i = 0; i < n; i++)
	{
		resc[i] = secant_method(x1[i],x2[i], epsilon);
		printf("初值 = %.14lf, %.14f, 根 = %.14lf, 迭代步数 = %d\n", x1[i], x2[i], resc[i].x0, resc[i].steps);
	}
	printf("\n");

	//验证牛顿法二阶收敛
	printf("验证牛顿法二阶收敛\n");
	for (i = 0; i < m; i++)
	{
		printf("初始值 = %.14lf\n", x0[i]);
		for (j = 1; j <= resn[i].steps; j++)
		{
			printf("(x* - x%d)/(x* - x%d)^2 = %.14lf\n", j, j - 1, fabs(resn[i].xk[j] - resn[i].x0) / (resn[i].xk[j - 1] - resn[i].x0) / (resn[i].xk[j - 1] - resn[i].x0));
		}
		printf("\n");
	}

	printf("判断前两个初始值是否为三阶收敛\n");
	for (i = 0; i < 2; i++)
	{
		printf("初始值 = %.14lf\n", x0[i]);
		for (j = 1; j <= resn[i].steps; j++)
		{
			printf("(x* - x%d)/(x* - x%d)^3 = %.14lf\n", j, j - 1, fabs(resn[i].xk[j] - resn[i].x0) / (resn[i].xk[j - 1] - resn[i].x0) / (resn[i].xk[j - 1] - resn[i].x0) / (resn[i].xk[j - 1] - resn[i].x0));
		}
		printf("\n");
	}

	printf("可以看出在有些特殊的初始值下，牛顿迭代可以达到三阶收敛\n");

	//验证弦截法二阶收敛
	printf("验证弦截法法(sqrt(5) + 1) / 2阶收敛\n");
	for (i = 0; i < n; i++)
	{
		printf("初始值 = %.14lf, %.14lf\n", x1[i], x2[i]);
		for (j = 2; j <= resc[i].steps + 1; j++)
		{
			printf("(x* - x%d)/(x* - x%d)^((sqrt(5) + 1) / 2) = %.14lf\n", j, j - 1, fabs(resc[i].xk[j] - resc[i].x0) / pow(fabs(resc[i].xk[j - 1] - resc[i].x0), (sqrt(5) + 1) / 2));
		}
		printf("\n");
	}

	//输出计算列表
	printf("牛顿法\n");
	for (i = 0; i < m; i++)
	{
		printf("初始值 = %.14lf\n", x0[i]);
		for (j = 0; j <= resn[i].steps; j++)
		{
			printf("k = %d\t xk = %+.14lf\t f(xk) = %+.14lf7\t error = %.14lf\n", j, resn[i].xk[j], resn[i].fk[j], resn[i].error[j]);
		}
		printf("x0 = %.14lf\n\n", resn[i].x0);
	}

	printf("弦截法\n");
	for (i = 0; i < n; i++)
	{
		printf("初始值 = %.14lf, %.14lf\n", x1[i], x2[i]);
		for (j = 0; j <= resc[i].steps + 1; j++)
		{
			printf("k = %d\t xk = %+.14lf\t f(xk) = %+.14lf7\t error = %.14lf\n", j, resc[i].xk[j], resc[i].fk[j], resc[i].error[j]);
		}
		printf("x0 = %.14lf\n\n", resc[i].x0);
	}

}

double f(double x)
{
	return x * x * x / 3 - x;
}

double f1(double x)
{
	return x * x - 1;
}

//定义牛顿法
result newton_method(double x0, double epsilon)
{
	int j = 0, k = 0, counter = 2;
	result res = {};
	double* error = (double *)malloc(N * sizeof(double));
	double* x = (double *)malloc(N * sizeof(double));
	double* fk = (double *)malloc(N * sizeof(double));
	x[0] = x0;
	fk[0] = f(x0);
	error[0] = -1;
	do
	{
		//防止出现列表长度不够的情况，使用realloc重新分配内存
		if (j >= N - 1)
		{
			error = (double *)realloc(error, counter * N * sizeof(double));
			x = (double *)realloc(x, counter * N * sizeof(double));
			fk = (double *)realloc(fk, counter * N * sizeof(double));
			if (!error || !x || ! fk)
			{
				free(error);
				free(fk);
				free(x);
				printf("realloc申请内存失败\n");
				return res;
			}
			counter++;
			j = 0;
		}
		j++;
		k++;
		//牛顿法公式
		x[k] = x[k - 1] - f(x[k - 1]) / f1(x[k - 1]);
		fk[k] = f(x[k]);
		error[k] = fabs(x[k] - x[k - 1]);

	} while (error[k] >= epsilon || fabs(fk[k]) >= epsilon);//判断满足条件，退出循环
	res.steps = k;//牛顿法有k + 1个点， 最后一个点的下标为k
	res.xk = x;
	res.fk = fk;
	res.error = error;
	res.x0 = x[k];
	return res;
}

//定义弦截法
result secant_method(double x0, double x1, double epsilon)
{
	int j = 1, k = 1, counter = 2;
	result res = {};
	double* error = (double *)malloc(N * sizeof(double));
	double* x = (double *)malloc(N * sizeof(double));
	double* fk = (double *)malloc(N * sizeof(double));
	x[0] = x0;
	x[1] = x1;
	fk[0] = f(x0);
	fk[1] = f(x1);
	error[0] = error[1] = -1;
	do
	{
		//防止出现列表长度不够的情况，使用realloc重新分配内存
		if (j >= N - 1)
		{
			error = (double *)realloc(error, counter * N * sizeof(double));
			x = (double *)realloc(x, counter * N * sizeof(double));
			fk = (double *)realloc(fk, counter * N * sizeof(double));
			if (!error || !x || !fk)
			{
				free(error);
				free(fk);
				free(x);
				printf("realloc申请内存失败\n");
				return res;
			}
			counter++;
			j = 0;
		}
		j++;
		k++;
		// 弦截法递推公式
		x[k] = x[k - 1] - f(x[k - 1])* (x[k - 1] - x[k - 2]) / (f(x[k - 1]) - f(x[k - 2]));
		fk[k] = f(x[k]);
		error[k] = fabs(x[k] - x[k - 1]);

	} while (error[k] >= epsilon || fabs(fk[k]) >= epsilon);//判断满足条件，退出循环
	res.steps = k - 1;//弦截法有k+1个点，进行了k-1次迭代，最后一个点的下标为k
	res.xk = x;
	res.fk = fk;
	res.error = error;
	res.x0 = x[k];
	return res;
}