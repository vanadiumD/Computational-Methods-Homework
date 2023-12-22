#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double f(double x);
int main()
{
	// 初始化变量
	// n为二分次数，初始为三个差值点
	// 按定义第 2^i + 1个结点积分应该有 N - i个误差阶，第一次积分无法计算误差，没有误差阶
	// o_ik 满足如下排列， i为第i次积分， k为分点的倍数:
	// o_2,2, o_2,4, ..., o_2,2^(N - 3), o_2,2^(N - 2)
	// o_3,2, o_3,4 ...,  o_3,2^(N - 3)
	//     ...
	// o_(N - 1),2
	int N = 12;
	double a = 1, b = 5;
	double h, sum;
	int n, i, j;
	double * T = (double *)malloc(N * sizeof(double));
	double * S = (double *)malloc(N * sizeof(double));
	double * et = (double *)malloc(N * sizeof(double));
	double * es = (double *)malloc(N * sizeof(double));
	double * ot = (double *)malloc((N - 2) * (N - 1) / 2 * sizeof(double));
	double * os = (double *)malloc((N - 2) * (N - 1) / 2 * sizeof(double));

	// 计算复化梯形积分并计算误差
	h = (b - a) / 2;
	j = 2;
	T[0] = h * (f(a) / 2 + f((a + b) / 2) + f(b) / 2);
	et[0] = -1;//误差从T_2算起，-1占位
	for (n = 1; n < N; n++)
	{
		h /= 2;
		j = j << 1;
		sum = T[n - 1] / 2;
		for (i = 1; i < j; i += 2)
		{
			sum += h * f(a + i * h);
		}
		T[n] = sum;
		//计算误差
		et[n] = fabs(sum - T[n - 1]);
	}
	//计算误差阶
	// 第i行
	for (i = 1; i < N ; i++)
	{
		n = 1;
		// 第j列 (第i行前面已经有(N - 2 + N - i) * (i - 1) / 2个元素)
		for (j = 1; j < N  - i; j++)
		{
			n = n << 1;
			ot[((2 * N - 2 - i) * (i - 1)) / 2 + j - 1] = log(et[i] / et[i + j]) / log(n);
		}
	}

	//输出结果
	printf("复化梯形积分\n");
	printf("T_1 = %.14lf\t",  T[0]);
	printf("e_h = %.14lf\n\n", et[0]);
	for (i = 1; i < N ; i++)
	{
		printf("T_%d = %.14lf\t", i + 1, T[i]);
		printf("e_h = %.14lf\n", et[i]);
		n = 1;
		// 第j列 第i行前面已经有(N - 1 + N - 1 - i) * i / 2个元素
		for (j = 1; j < N  - i; j++)
		{
			n = n << 1;
			printf("o_%d = %.8lf\t", n, ot[((2 * N - 2 - i) * (i - 1)) / 2 + j - 1]);
		}
		printf("\n\n");
	}

	// 计算复化辛普森积分并计算误差
	h = (b - a) / 2;
	j = 2;
	S[0] = h * (f(a) + 4 * f((a + b) / 2) + f(b)) / 3;
	es[0] = -1;//误差从T_2算起，-1占位
	for (n = 1; n < N; n++)
	{
		sum = S[n - 1] / 2;
		// 将原来系数为4/3的奇数插值点变为系数为2/3的偶数插值点
		for (i = 1; i < j; i += 2)
		{
			sum -= f(a + i * h) * h / 3;
		}
		j = j << 1;
		h /= 2;
		// 加入新的系数为4/3的奇数插值点
		for (i = 1; i < j; i += 2)
		{
			sum += f(a + i * h) * 4 * h / 3;
		}
		S[n] = sum;
		//计算误差
		es[n] = fabs(sum - S[n - 1]);
	}
	//计算误差阶
	// 第i行
	for (i = 1; i < N; i++)
	{
		n = 1;
		// 第j列 (第i行前面已经有(N - 2 + N - i) * (i - 1) / 2个元素)
		for (j = 1; j < N - i; j++)
		{
			n = n << 1;
			os[((2 * N - 2 - i) * (i - 1)) / 2 + j - 1] = log(es[i] / es[i + j]) / log(n);
		}
	}

	//输出结果
	printf("复化辛普森积分\n");
	printf("S_1 = %.14lf\t", S[0]);
	printf("e_s = %.14lf\n\n", es[0]);
	for (i = 1; i < N; i++)
	{
		printf("S_%d = %.14lf\t", i + 1, S[i]);
		printf("e_s = %.14lf\n", es[i]);
		n = 1;
		// 第j列 第i行前面已经有(N - 1 + N - 1 - i) * i / 2个元素
		for (j = 1; j < N - i; j++)
		{
			n = n << 1;
			printf("o_%d = %.8lf\t", n, os[((2 * N - 2 - i) * (i - 1)) / 2 + j - 1]);
		}
		printf("\n\n");
	}

	//释放内存
	free(T);
	free(et);
	free(ot);
	free(S);
	free(es);
	free(os);
}

double f(double x)
{
	return sin(x);
}

