#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//限定矩阵范围
#define N 20
#define M 20
//自定义矩阵运算函数
//矩阵计算函数由本人于大一计算机程序设计A中大作业原创
void matrix_multiply(double *a, double *b, int m, int s, int n, double *res);//矩阵乘法
bool initmatrix(double *a, int m, int n);//初始化矩阵(储存结果)
bool identity_matrix(double *a, int n);//创建单位矩阵
void matrix_copy(double *a, int m, int n, double * res);//拷贝矩阵(避免覆盖)
bool submatrix(double *a, int m, int n, int row_start, int row_length, int column_start, int column_length, double *res);//提取子矩阵
int solve_liner_equation(double *a, int m, int n, double * res);//利用高斯消元求解方程组
int row_reduce(double *a, int m, int n, double *res, char simp, int model);//行初等变换化为上三角矩阵
int least_square_method(double *a, int m, int n, double *res);//最小二乘法
void output_matrix(double *a, int m, int n);//输出矩阵
void transform(double *a, int m, int n, double *res);//转置矩阵

int main()
{
	//输入列表
	int i, m, n;
	double x_list[] = { 0.25, 0.5, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50 };
	double y_list[] = { 1.284, 1.648, 2.117, 2.718, 3.427, 2.798, 3.534, 4.456, 5.465, 5.894 };
	double error, a, b;
	//计算矩阵A的元素(sin(x), cos(x), y)
	m = sizeof(x_list) / sizeof(x_list[0]);
	n = 2;
	double * A = (double*)malloc(m * (n + 1) * sizeof(double));
	double * res = (double*)malloc(n * sizeof(double));
	for (i = 0; i < m; i++)
	{
		A[i * (n + 1)] = sin(x_list[i]);
	}
	for (i = 0; i < m; i++)
	{
		A[i * (n + 1) + 1] = cos(x_list[i]);
	}
	for (i = 0; i < m; i++)
	{
		A[i * (n + 1) + 2] = y_list[i];
	}
	//利用最小二乘法函数计算
	least_square_method(A, m, n + 1, res);
	a = res[0], b = res[1];
	//计算均方误差：
	error = 0;
	for (i = 0; i < m; i++)
	{
		error += pow((A[i * (n + 1)] * a + A[i * (n + 1) + 1] * b - A[i * (n + 1) + 2]), 2);
	}
	error /= m;
	printf("a = %lf, b = %lf, 均方误差 = %lf", a, b, error);
}
//最小二乘法：a是(A,b)型增广矩阵，通过提取子矩阵A和b，计算A^T*(A,b),可得(A^T A, A^T,B)，再对增广矩阵进行变换
//化为上三角矩阵，再利用高斯消元即可求出最终的解
int least_square_method(double *a, int m, int n, double *res)
{
	//创建临时矩阵储存A转置
	double *at, *temp;
	at = (double*)malloc(m * (n - 1) * sizeof(double));
	temp = (double*)malloc(m * n * sizeof(double));
	int r;
	//提取增广矩阵a的前n - 1列作为A
	submatrix(a, m, n, 0, m, 0, n - 1, at);
	transform(at, m, n - 1, at);
	//计算A^T a = A^T (A, b)
	matrix_multiply(at, a, n - 1, m, n, temp);
	//解增广矩阵对应的方程组,r为增广矩阵的秩
	r = solve_liner_equation(temp, n - 1, n, res);
	free(at);
	free(temp);
	return r;
}

void matrix_multiply(double *a, double *b, int m, int s, int n, double *res)
{
	int i, j, k;
	double *temp1, *temp2;
	temp1 = (double *)malloc(m * s * sizeof(double));
	temp2 = (double *)malloc(s * n * sizeof(double));
	matrix_copy(a, m, s, temp1);
	matrix_copy(b, s, n, temp2);
	if (!initmatrix(res, m, n))
	{
		printf("初始化结果矩阵失败");
		return;
	}

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < s; k++)
			{
				res[i * n + j] += temp1[i * s + k] * temp2[k * n + j];
			}
		}
	}
	free(temp1);
	free(temp2);
}

bool initmatrix(double *a, int m, int n)
{
	int i, j;

	if (m > N || n > N)
		return false;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i * n + j] = 0;
		}
	}
	return true;
}

void matrix_copy(double *a, int m, int n, double * res)
{
	int i, j;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			res[i * n + j] = a[i * n + j];
		}
	}
}

void output_matrix(double *a, int m, int n)
{
	int i, j;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			printf("% lf\t", a[i * n + j]);
		}
		printf("\n");
	}
}

bool identity_matrix(double *a, int n)
{
	int i, j;

	if (n > M || n > N)
		return false;

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			a[i * n + j] = 0;
		}
	}

	for (i = 0; i < n; i++)
		a[i * n + i] = 1;
	return true;
}

int solve_liner_equation(double *a, int m, int n, double * res)
{
	int i, j, k, dim, rank, pos[N];
	double judge = 0;

	//创建temp矩阵, 将a化为简化阶梯型矩阵后存在temp中
	double *temp;
	temp = (double*)malloc(m * n * sizeof(double));
	if ((rank = row_reduce(a, m, n, temp, 's', 0)) == 0)
	{
		free(temp);
		if (identity_matrix(res, n) == false)
			return -1;
		else
			return n;
	}

	//无解
	for (i = 0; i < n - 1; i++)
		judge += fabs(temp[(rank - 1) * n + i]);
	if (judge < 1e-20)
	{
		free(temp);
		return 0;
	}

	//求特解, 解空间向量维数为n - 1
	//res大小为(n - 1) * (n - rank);
	if (initmatrix(res, n - 1, n - rank) == false)
	{
		free(temp);
		return -1;
	}
	//矩阵列满秩, 只有唯一解
	if (rank == n - 1)
	{
		for (i = 0; i < n - 1; i++)
			res[i] = temp[n * i + n - 1];
		free(temp);
		return 1;
	}

	//有多个解:
	//求特解
	dim = n - rank;
	k = 0;
	for (i = 0; i < N; i++)
		pos[i] = 0;
	for (i = 0, j = 0; i < m && j < n - 1; i++, j++)
	{
		while (j < n - 1 && fabs(temp[i * n + j]) < 1e-20)
		{
			res[j * dim] = 0;
			pos[k++] = j++;
		}
		res[j * dim] = temp[i * n + n - 1];
	}
	//将行满秩最后一行主元以后的自由变量记录位置
	while (j < n - 1)
		pos[k++] = j++;

	//求通解, 齐此方程解空间的维数为dim - 1;
	for (k = 1; k < dim; k++)
	{
		for (i = 0, j = 0; i < m && j < n - 1; i++, j++)
		{
			while (fabs(temp[i * n + j]) < 1e-20)
			{
				res[j * dim + k] = 0;
				j++;
			}
			res[j * dim + k] = -temp[i * n + pos[k - 1]];
		}
		res[pos[k - 1] * dim + k] = 1;

	}
	free(temp);
	return dim;
}

static int factor;
int row_reduce(double *a, int m, int n, double *res, char simp, int model)
{
	matrix_copy(a, m, n, res);

	int i, j, k, s, t, flag, rank = 0, step = 1;
	double temp;
	factor = 1;
	for (i = 0, j = 0; i < m && j < n; j++)
	{
		//找j列非零元
		k = i;
		while (fabs(res[k * n + j]) < 1e-20 && ++k < m);
		if (k == m)
			continue;

		rank++;//存在一个零元, 矩阵的秩加1

		//找到最后一行的非零元后退出消元
		if (i == m - 1 && simp != 's')
			break;

		//将有非零元的一行交换到i行
		if (i != k)
		{
			for (s = j; s < n; s++)
			{
				temp = res[k * n + s];
				res[k * n + s] = res[i * n + s];
				res[i * n + s] = temp;
			}
			if (model == 1)
			{
				printf("step%d(S):\n", step++);
				output_matrix(res, m, n);
			}
			factor *= -1;//每次交换全局变量factor变号;
		}

		//利用主元将j列其他元素消为零
		flag = 0;
		if (simp == 's')
			s = 0;
		else s = i + 1;
		for (; s < m; s++)
		{
			if (s == i || fabs(temp = res[s * n + j] / res[i * n + j]) < 1e-20)
				continue;//保存主元乘的系数系数
			for (t = j; t < n; t++)
			{
				flag = 1;
				res[s * n + t] -= temp * res[i * n + t];
			}
		}

		if (model == 1 && flag == 1)
		{
			printf("step%d(T):\n", step++);
			output_matrix(res, m, n);
		}

		if (simp == 's')
		{
			temp = res[i * n + j];
			for (s = j; s < n; s++)
				res[i * n + s] /= temp;
		}

		i++;
	}
	return rank;
}

bool submatrix(double *a, int m, int n, int row_start, int row_length, int column_start, int column_length, double *res)
{

	if (m < row_start + row_length || n < column_start + column_length)
		return false;
	int i, j;
	a += row_start * n + column_start;

	for (i = 0; i < row_length; i++)
	{
		for (j = 0; j < column_length; j++)
		{
			res[i * column_length + j] = a[i * n + j];
		}
	}
	return true;
}

void transform(double *a, int m, int n, double *res)
{

	int i, j;
	double *temp;
	temp = (double *)malloc(m * n * sizeof(double));
	matrix_copy(a, m, n, temp);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			res[j * m + i] = temp[i * n + j];
		}
	}
	free(temp);
}

