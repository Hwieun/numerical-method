#include "stdafx.h"
#include<math.h>
#include<stdlib.h>

#define TOL 0.000001

typedef struct {
	int test;
}Thread;
typedef struct {
	Thread* head;
}List;

double p_data[10];
double F(double t, double y) {
	return 10 * exp(5 * t)*(y - t);
}

double f(double t, double y) {
	return 5 * exp(5 * t)*(y - t)*(y- t) + 1;
}
double abs(double result)
{
	if (result < 0)
		result = result*(-1);
	return result;
}
double exact(double t) {
	return t - exp((-5)*t);
}
double Trapeziodal(double y, double h, double t) {
	double ykNext = 0;
	double yk = y;
	char M = 0;
	ykNext = yk - (yk - y - h*0.5*(f(t, y) + f(t + h, yk))) / (1 - h*0.5*F(t + h, yk));

	while ((abs(ykNext - yk) > TOL)&&(M<10)) {
		yk = ykNext;
		ykNext = yk - (yk - y - h*0.5*(f(t, y) + f(t + h, yk))) / (1 - h*0.5*F(t + h, yk));
		
		M++;
	}
	return y + h*0.5*(f(t + h, ykNext) + f(t, y));
}
void cal(double h) {
	double y = -1; double t0 = 0;

	for (int i = 0; i <6; i++) {
		p_data[i] = y;

		y = Trapeziodal(y, h, t0);
		t0 += h;
	}
}
double Error(double t, int argu) {
	double result = exact(t) - p_data[argu];
	return abs(result);
}
int main() {
	int count = 0;

	printf("       Trapezoidal h=0.2\n");
	printf("ti              wi                | y(ti)-wi |\n");
	printf("----------------------------------------------\n");

	cal(0.2);

	for (int i = 0; i < 6; i++) {
		printf("%.1lf       %.7lf                %.8lf\n", i*0.2, p_data[i], Error(0.2*i, i));
	}
	cal(0.25);
	printf("h=0.25\n");
	for (int i = 0; i < 5; i++) {
		printf("%.2lf       %.7lf                %.8lf\n", i*0.25, p_data[i], Error(0.25*i,i));
	}
	
	List* list = (List*)malloc(sizeof(List));
	Thread* TCB = list->head;
	
	system("pause");

	return 0;
}