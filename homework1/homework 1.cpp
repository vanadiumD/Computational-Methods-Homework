#include<stdio.h>
#include<math.h>
#define PI acos(-1)
int main() {
	double x[] = { 0.0, 0.5, 1.0, sqrt(2), 10.0, 100.0, 300.0 };
	double a, b, ans, res;
	int i;
	for (i = 0; i < 7; i++)
	{
		ans = 0;
		res = 0;
		a = 1;
		b = 1 + x[i];
		//利用不等式 1/k(k + x) <= 1/k^2, 以及级数sum_0^infty 1/k^2 = pi^2 / 6 
		//可以得到余项的不等式 res(1/k(k + x)) <= res(1/k^2) = pi^2 / 6  - res
		//其中，res 为1/k^2 的前n项和
		while (fabs(PI * PI / 6 - res) >= 1e-7) 
		{
			ans += 1 / (a * b);
			res += 1 / (a * a);
			a++;
			b++;
		}
		printf("x=%.2lf,y=%.15lf\n", x[i], ans);
	}
	return 0;
}