#include "stdafx.h"
#include<math.h>
#include<stdlib.h>
#define TOL 0.000001
double p_data[2][10];
//int count = 0;
double f(double t, double y) {
	return 5 * exp(5 * t)*(y - t)*(y - t) + 1;
}
double abs(double result)
{
	if (result < 0)
		result = result*(-1);
	return result;
}
double RungeKutta(double y, double h, double t) {
	double k1 = h*f(t, y);
	double k2 = h*f(t + h*0.5, y + k1*0.5);
	double k3 = h*f(t + h*0.5, y + k2*0.5);
	double k4 = h*f(t+h, y + k3);
	return y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}
double exact(double t) {
	return t - exp((-5)*t);
}
/*double F(double t, double y) {
return 10 * exp(5 * t)*(y - t);
}
double Trapeziodal(double y, double h, double t) {
double yNext = y - (y - h*0.5*f(t + h, y)) / (1 - 0.5*h*F(t + h, y));
if (abs(yNext - y) < TOL)
return yNext;
}*/
void cal(double fp(double ,double, double), double h, int argu) {
	double y = -1; double t0 = 0;
	int N = 1 / h; //h=0.2, N=5
	int flag;
	if (fp == RungeKutta) flag = 1;
	else flag = 2;
	for (int i = 0; i <= N*flag; i++) {
	// if (i%2== 0)
		p_data[argu][i/flag] = y;
		y = fp(y, h, t0);
		t0 += h;
	}
}
double Error(double t, int argu1, int argu2) {
	double result = exact(t) - p_data[argu1][argu2];
	return abs(result);
}
int main() {
	int count = 0;
	printf(" Runge-Kutta Order4 h=0.2\n");
	printf("ti wi | y(ti)-wi |\n");
	printf("----------------------------------------------\n");
	cal(RungeKutta, 0.2, 0);
	//cal(Trapeziodal, 0.2, 1);
	for (int i = 0; i < 6; i++) {
		printf("%.1lf %.7lf %.8lf\n", i*0.2, p_data[0][i],
		Error(0.2*i, 0, i));
	}
	cal(RungeKutta, 0.25, 0);
	//cal(Trapeziodal, 0.25, 1);
	printf("h=0.25\n");
	for (int i = 0; i < 4; i++) {
		printf("%.2lf %.7lf %.8lf\n", i*0.25, p_data[0][i], Error(0.25*i, 0, i));
	}
	printf("%.2lf Overflow %.8lf\n", 1, Error(1, 0, 4));
	return 0;
}