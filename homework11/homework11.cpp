#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define N 2//定义问题维度为2

double h = 1e-5;//定义计算梯度和势能矩阵步长为1e-5
double f(double *x)//定义目标函数
{
	return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
}
double grad_f(double *x, double *grad_x);//计算目标函数梯度
void hessian_f(double *x, double *G);//计算hessian矩阵
void inverse(double *a, double *b);//利用高斯-若尔当法求逆矩阵
void multiply(double *a, double *x, double *b);//计算矩阵乘法
void steepest_descent_method(double *x, double error);//最速下降法
void newton_method(double *x, double error);//牛顿法
void backtracing_method(double *x, double *p, double h, double lambda, double *a, double *b); //前进法
void golden_section_search(double *x, double *p, double a, double b, double error);//黄金分割一维搜索法
void output(int i, double *x);//输出结果函数

int main()
{
	double x1[2] = { 0, 0 }, x2[2] = { 0, 0 };
	double error = 1e-4;
	steepest_descent_method(x1, error);
	newton_method(x2, error);
	//double a[2] = { 5, 6}, b[4], c[4];
	//hessian_f(a, b);
	//inverse(b, c);
	//printf("a");
}

void steepest_descent_method(double *x, double error)
{
	double *p = (double*)malloc(N * sizeof(double));
	double e, a, b;
	int k = 0;
	do
	{
		e = grad_f(x, p);
		backtracing_method(x, p, 0.3, 0, &a, &b);
		golden_section_search(x, p, a, b, error);
		k++;
		output(k, x);
	} while (e > error);
}

void newton_method(double *x, double error)
{
	int i, k;
	double *a, *b, *c;
	double e;
	a = (double*)malloc(N * N * sizeof(double));
	b = (double*)malloc(N * N * sizeof(double));
	c = (double*)malloc(N * sizeof(double));
	k = 0;
	do
	{
		hessian_f(x, b);
		inverse(b, a);
		e = grad_f(x, b);
		multiply(a, b, c);
		for (i = 0; i < N; i++)
		{
			x[i] -= c[i];
		}
		k++;
		output(k, x);
	} while (e > error);
}

void backtracing_method(double *x, double *p, double h, double lambda, double *a, double *b)
{
	int i, k;
	double fa, fb, *x1, lambda1;
	x1 = (double *)malloc(N * sizeof(double));
	for (i = 0; i < N; i++)
	{
		x1[i] = x[i] + lambda * p[i];
	}
	fa = f(x1);
	k = 0;
	while(1)
	{
		lambda1 = lambda + h;
		for (i = 0; i < N; i++)
		{
			x1[i] = x[i] + lambda1 * p[i];
		}
		fb = f(x1);
		if (0 == k && fa < fb)
		{
			h = -h;
			lambda = lambda1;
			k++;
			continue;
		}
		if (fa < fb && 0 != k)
		{
			if (lambda < lambda1)
			{
				*a = lambda;
				*b = lambda1;
			}
			else
			{
				*a = lambda1;
				*b = lambda;
			}
			free(x1);
			return;
		}
		h *= 2;
		lambda = lambda1;
		k++;
	}
}

void golden_section_search(double *x, double *p, double a, double b, double error)
{
	int i;
	double alpha, beta, tau, e;
	double *x1, *x2;
	x1 = (double *)malloc(N * sizeof(double));
	x2 = (double *)malloc(N * sizeof(double));
	tau = (sqrt(5) - 1) / 2;
	alpha = b - tau * (b - a);
	beta = a + tau * (b - a);
	do
	{
		for (i = 0; i < N; i++)
		{
			x1[i] = x[i] + alpha * p[i];
			x2[i] = x[i] + beta * p[i];
		}
		if (f(x1) > f(x2))
		{
			a = alpha;
			alpha = beta;
			beta = a + tau * (b - a);
		}
		else
		{
			b = beta;
			beta = alpha;
			alpha = b - tau * (b - a);
		}
		e = b - a;

	} while (e > error);
	for (i = 0; i < N; i++)
	{
		x[i] = x[i] + (a + b) / 2 * p[i];
	}
	free(x1);
	free(x2);
}

double grad_f(double *x, double *grad_x)
{
	int i;
	double model = 0;
	for (i = 0; i < N; i++)
	{
		//利用中心差商计算梯度
		x[i] += h;
		grad_x[i] = f(x);
		x[i] -= 2 * h;
		grad_x[i] -= f(x);
		x[i] += h;
		grad_x[i] /= 2 * h;
		model += grad_x[i] * grad_x[i];
	}
	return sqrt(model);
}

void hessian_f(double *x, double *G)
{
	int i, j;
	double ans, fx;
	for (i = 0; i < N; i++)
	{
		//利用五点差分法计算对角元
		x[i] += 2 * h;
		ans = -f(x);
		x[i] -= h;
		ans += 16 * f(x);
		x[i] -= h;
		ans -= 30 * f(x);
		x[i] -= h;
		ans += 16 * f(x);
		x[i] -= h;
		ans -= f(x);
		x[i] += 2 * h;
		G[i * N + i] = ans / (12 * h * h);

		//利用中心差商计算非对角元
		for (j = i + 1; j < N; j++)
		{
			x[i] += h;
			x[j] += h;
			ans = f(x);
			x[j] -= 2 * h;
			ans -= f(x);
			x[i] -= 2 * h;
			ans += f(x);
			x[j] += 2 * h;
			ans -= f(x);
			G[i * N + j] = ans / (4 * h * h);
			G[j * N + i] = G[i * N + j];
			x[i] += h;
			x[j] -= h;
		}
	}
}

void inverse(double * a, double *b)
{
	int i, j, k, s;
	double max, temp;
	for (i = 0; i < N; i++)
	{
		b[i * N + i] = 1;
		for (j = i + 1; j < N; j++)
		{
			b[i * N + j] = 0;
			b[j * N + i] = 0;
		}
	}

	for (i = 0; i < N - 1; i++)
	{
		s = i;
		max = fabs(a[i * N + i]);
		for (k = i + 1; k < N; k++)
		{
			if (max < fabs(a[k * N + i]))
			{
				max = fabs(a[k * N + i]);
				s = k;
			}
		}

		if (s != i)
		{
			for (k = i; k < N; k++)
			{
				temp = a[i * N + k];
				a[i * N + k] = a[s * N + k];
				a[s * N + k] = temp;
			}
			for (k = 0; k < N; k++)
			{
				temp = b[i * N + k];
				b[i * N + k] = b[s * N + k];
				b[s * N + k] = temp;
			}
		}

		for (j = i + 1; j < N; j++)
		{
			a[j * N + i] /= a[i * N + i];
			for (k = i + 1; k < N; k++)
			{
				a[j * N + k] -= a[i * N + k] * a[j * N + i];
			}
			for (k = 0; k < N; k++)
			{
				b[j * N + k] -= b[i * N + k] * a[j * N + i];
			}
		}
	}

	for (j = 0; j < N; j++)
	{
		for (i = N - 1; i >= 0; i--)
		{
			for (k = i + 1; k < N; k++)
			{
				b[i * N + j] -= b[k * N + j] * a[i * N + k];
			}
			b[i * N + j] /= a[i * N + i];
		}
	}

}

void multiply(double *a, double *x, double *b)
{
	int i, j;
	for (i = 0; i < N; i++)
	{
		b[i] = 0;
		for (j = 0; j < N; j++)
		{
			b[i] += a[i * N + j] * x[j];
		}
	}
}

void output(int i, double *x)
{
	printf("第%d次迭代，f(x_%d) = %lf, x_%d = %lf, y_%d = %lf\n", i, i, f(x), i, x[0], i, x[1]);
}