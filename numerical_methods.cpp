// hw8.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <math.h>
double g_data[3][6];
double t0, y, h;
typedef double(*FP)(double, double, double);
double y_prime(double t, double y) {
	return y - t*t + 1;
}
double Euler(double y, double h, double t) {
	return y + h*y_prime(t, y);
}
double M_Euler(double y, double h, double t) {
	return y + (h / 2.0)*(y_prime(t, y) + y_prime(t+h, y + h*y_prime(t,y)));
}
double RungeKutta(double y, double h, double t) {
	double k1, k2, k3, k4;
	k1 = h*y_prime(t, y);
	k2 = h*y_prime(t + h / 2.0, y + k1 / 2.0);
	k3 = h*y_prime(t + h / 2.0, y + k2 / 2.0);
	k4 = h*y_prime(t + h, y + k3);
	return y + (1 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4);
}
double Exact(double t) {
	return (t + 1)*(t + 1) - 0.5*exp(t);
}
void init(double *p_t, double *p_y) {
	*p_t = 0;
	*p_y = 0.5;
}
void numerical(FP fp, double *p_data, double a_h) {
	init(&t0, &y);
	h = a_h;
	int filter = 1;
	while (a_h*filter != 0.1) filter++;
	for (int i = 0; i <= filter * 5; i++) {
		if (i%filter == 0) {
			p_data[i/filter] = y;
		}
		y = fp(y, h, t0);
		t0 += h;
	}
}
void print() {
	printf("ti Exact Euler Modified Euler Runge-Kutta Order 4\n");
	double t = 0;
	for (int i = 0; i <= 5; i++) {
		printf("%.1lf %.7lf %.7lf %.7lf %.7lf\n", t, Exact(t), g_data[0][i], g_data[1][i], g_data[2][i]);
		t += 0.1;
	}
}
void calculate() {
	FP table[3] = { Euler, M_Euler, RungeKutta };
	h = 0.025;
	for (int i = 0; i < 3; i++) {
		numerical(table[i], g_data[i], h);
		h += h;
	}
	print();
}
int main()
{
	calculate();
	return 0;
}