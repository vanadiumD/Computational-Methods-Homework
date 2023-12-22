#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define M 4

//储存结果的结构体
typedef struct RESULT
{
	double *x = NULL;
	double *y = NULL;
	double h = 0;
	int N = 0;
}result;

//一阶微分方程dy/dx = f(x, y)
//f(x, y)表达式
double ode(double x, double y)
{
	return -x * x * y * y;
}

//精确解函数
double f(double x)
{
	return 3 / (1 + x * x * x);
}

result runge_kutta4(double ode(double, double), double x0, double x1, double y0, double h);//四阶龙格库塔法
result implicit_adams4(double ode(double, double), double x0, double x1, double y0, double h);//四阶隐式adams方法
void print(result res, int N);//输出结果结构体中的x,y，用于debug

int main()
{
	int i, j;
	double error1[M], error2[M];
	result res;
	double x0 = 0, x1 = 1.5, y0 = 3, ans = f(x1), h;

	//利用两种方法求解
	for (i = 0; i < M; i++)
	{
		h = 0.1 / pow(2, i);
		res = runge_kutta4(ode, x0, x1, y0, h);
		//print(res, res.N + 1);
		error1[i] = fabs(res.y[res.N] - ans);
		res = implicit_adams4(ode, x0, x1, y0, h);
		//print(res, res.N + 1);
		error2[i] = fabs(res.y[res.N] - ans);
	}
	//计算误差和误差阶
	printf("四阶Runge_Kutta公式的误差和误差阶\n");
	for (i = 0; i < M - 1; i++)
	{
		h = 0.1 / pow(2, i);
		printf("h = %lf, err = %.12lf", h, error1[i]);
		for (j = i + 1; j < M; j++)
		{
			printf(", o%d = %.12lf", int(pow(2, j - i)), log(error1[i] / error1[j]) / log(2) / (j - i));
		}
		printf("\n");
	}
	printf("\n");

	printf("四阶隐式Adams公式的误差和误差阶\n");
	for (i = 0; i < M - 1; i++)
	{
		h = 0.1 / pow(2, i);
		printf("h = %lf, err = %.12lf", h, error2[i]);
		for (j = i + 1; j < M; j++)
		{
			printf(", o%d = %.12lf", int(pow(2, j - i)), log(error2[i] / error2[j]) / log(2) / (j - i));
		}
		printf("\n");
	}
}

result runge_kutta4(double ode(double, double), double x0, double x1, double y0, double h)
{
	int N = int(ceil((x1 - x0) / h));
	h = (x1 - x0) / N;

	int i;
	double *x, *y, k1, k2, k3, k4;
	result res;

	x = (double*)malloc((N + 1) * sizeof(double));
	y = (double*)malloc((N + 1) * sizeof(double));
	x[0] = x0;
	y[0] = y0;

	for (i = 0; i < N; i++)
	{
		x[i + 1] = x[i] + h;

		k1 = ode(x[i], y[i]);
		k2 = ode(x[i] + h / 2, y[i] + h * k1 / 2);
		k3 = ode(x[i] + h / 2, y[i] + h * k2 / 2);
		k4 = ode(x[i + 1], y[i] + h * k3);

		y[i + 1] = y[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;

	}
	res.y = y;
	res.x = x;
	res.h = h;
	res.N = N;
	return res;
}

result implicit_adams4(double ode(double, double), double x0, double x1, double y0, double h)
{
	int N = int(ceil((x1 - x0) / h));
	h = (x1 - x0) / N;
	double *x, *y, error, k1, k2, k3, k4, temp;
	result res;

	if (N < 3)
	{
		printf("隐式四阶adam方法最少需要3个结点\n");
		return res;
	}
	int i;

	x = (double*)malloc((N + 1) * sizeof(double));
	y = (double*)malloc((N + 1) * sizeof(double));
	x[0] = x0;
	y[0] = y0;

	//先用四阶龙格-库塔法计算前两个点
	for (i = 0; i < 3; i++)
	{
		x[i + 1] = x[i] + h;

		k1 = ode(x[i], y[i]);
		k2 = ode(x[i] + h / 2, y[i] + h * k1 / 2);
		k3 = ode(x[i] + h / 2, y[i] + h * k2 / 2);
		k4 = ode(x[i + 1], y[i] + h * k3);

		y[i + 1] = y[i] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
	}

	for (i = 3; i < N; i++)
	{
		x[i + 1] = x[i] + h;
		//计算显示adams公式
		temp = y[i] + (55 * ode(x[i], y[i]) - 59 * ode(x[i - 1], y[i - 1]) + 37 * ode(x[i - 2], y[i - 2]) - 9 * ode(x[i - 3], y[i - 3]))* h / 24;
		//计算隐式adams公式
		y[i + 1] = y[i] + (9 * ode(x[i + 1], temp) + 19 * ode(x[i], y[i]) - 5 * ode(x[i - 1], y[i - 1]) + ode(x[i - 2], y[i - 2])) * h / 24;
	}
	res.y = y;
	res.x = x;
	res.h = h;
	res.N = N;
	return res;
}

void print(result res, int N)
{
	int i;
	printf("n\t\tx\t\ty\n");
	for (i = 0; i < N; i++)
	{
		printf("%d\t%lf\t%lf\n", i, res.x[i], res.y[i]);
	}
}