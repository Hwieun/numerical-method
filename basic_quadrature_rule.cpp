#include "stdafx.h"
#include<stdlib.h>
#include <math.h>

double f(double x) {
	if (x == 1.8)
		return 3.12014;
	else if(x ==2.0)
		return 4.42569;
	else if (x == 2.2)
		return 6.04241;
	else if (x == 2.4)
		return 8.03014;
	else if (x == 2.6)
		return 10.46675;
}
double Midpoint(double a, double b) {
	double x0 = (a + b) / 2;
	return (b - a)*f(x0);
}
double Trapezoidal(double a, double b) {
	double result = (b - a)*(f(a) + f(b)) / 2;
	return result;
}
double Simson(double x0, double x2) {
	double x1 = (x0 + x2) / 2;
	return (x2 - x0)*(f(x0) + 4 * f(x1) + f(x2)) / 6;
}


int main() {

	printf("Midpoint Rule = %lf\n\n", Midpoint(1.8, 2.6));
	printf("Trapezoidal Rule = %lf\n\n", Trapezoidal(1.8, 2.6));
	printf("Simson's Rule = %lf\n", Simson(1.8, 2.6));

	getchar();

	return 0;
}
