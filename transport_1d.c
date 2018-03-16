#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

double A[100][100], B[100], X[100];
double dt = 0.1;
double dx = 1;
double v[100];
double lambda[100][3], delta[100], P1, P2, Q1, Q2, R1[100], R2[100], alpha[100][2];

double max3(double a, double b, double c)
{
	if ((a > b) && (a > c))
		return a;
	else if (b > c)
		return b;
	else
		return c;
}

double min2(double a, double b)
{
	if (a < b)
		return a;
	else
		return b;
}

double max2(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
}

double alpha_func(int j, int k)
{
	int a;
	if (j < k)
		a = j;
	else
		a = k;
	if (X[k] >= X[j])
		return min2(R1[k] * delta[a], lambda[j][1 + k - j] + delta[a]);
	else
		return min2(R2[k] * delta[a], lambda[j][1 + k - j] + delta[a]);
}

int main(void)
{
	int i, t, j;
	double end_T = 10;
	FILE *f;
	char file_name[20];
	double C[100], D[100];
	memset(X, 0, 100 * sizeof(double));
	for (i = 30; i < 40; i++)
		X[i] = 1;
	for (i = 0; i < 100; i++)
		v[i] = 1;
	printf("Init done\n");
	for (t = 0; t < end_T / dt; t ++) {
		printf("Time step %d\n", t);
		sprintf(file_name, "tmp/transport%d.dat", t);
		f = fopen(file_name, "w");
		for (i = 0; i < 100; i++)
			fprintf(f, "%10.10lf\n", X[i]);
		fclose(f);
		memset(A, 0, 10000 * sizeof(double));
		memset(B, 0, 100 * sizeof(double));
		/*	
		for (i = 0; i < 100; i++) {
			A[i][i] += 1.0 / dt;
			B[i] += X[i] / dt;
			if (i != 0) {
				A[i][i] += 1.0 / dx;
				A[i][i-1] += -1.0 / dx;
			}
		}
		for (i = 0; i < 100; i++) {
			if (i == 0)
				X[i] = B[i] / A[i][i];
			else
				X[i] = (B[i] - A[i][i-1] * X[i-1]) / A[i][i];
		}
		*/
		
		lambda[0][2] = - v[0 + 1] / (2 * dx);
		for (i = 1; i < 99; i++) {
			lambda[i][0] = v[i - 1] / (2 * dx);
			lambda[i][2] = - v[i + 1] / (2 * dx);
			lambda[i][1] = lambda[i][0] + lambda[i][2];
		}
		lambda[99][0] = v[99 - 1] / (2 * dx);
		for (i = 0; i < 99; i++) {
			delta[i] = max3(0, lambda[i][2], lambda[i + 1][0]);
		}
		for (i = 1; i < 99; i++) {
			Q1 = max2(0, lambda[i][0]) * max2(0, X[i - 1] - X[i]) + max2(0, lambda[i][2]) * max2(0, X[i +1] - X[i]);
			Q2 = max2(0, lambda[i][0]) * max2(0, X[i - 1] - X[i]) + min2(0, lambda[i][2]) * max2(0, X[i +1] - X[i]);
			P1 = min2(0, lambda[i][0]) * max2(0, X[i - 1] - X[i]) + min2(0, lambda[i][2]) * max2(0, X[i +1] - X[i]);
			P2 = min2(0, lambda[i][0]) * max2(0, X[i - 1] - X[i]) + max2(0, lambda[i][2]) * max2(0, X[i +1] - X[i]);
			if (P1 == 0)
				Q1 = 0;
			else
				Q1 /= P1;
			if (P2 == 0)
				Q2 = 0;
			else
				Q2 /= P2;
			R1[i] = max2(0, min2(1, Q1));
			R2[i] = max2(0, min2(1, Q2));
		}
		for (i = 1; i < 99; i++) {
			if (lambda[i][0] <= lambda[i - 1][2])
				alpha[i][0] = alpha_func(i, i - 1);
			else
				alpha[i][0] = alpha_func(i - 1, i);
			if (lambda[i][2] <= lambda[i + 1][0])
				alpha[i][1] = alpha_func(i, i + 1);
			else
				alpha[i][1] = alpha_func(i + 1, i);
		}
		for (i = 1; i < 99; i++) {
			lambda[i][0] += delta[i - 1];// - alpha[i][0];
			lambda[i][2] += delta[i];// - alpha[i][1];
		}
		A[0][0] += 1.0 / dt;
		B[0] += X[0] / dt;
		for (i = 1; i < 99; i++) {
			A[i][i] += 1.0 / dt;
			B[i] += X[i] / dt;
			A[i][i - 1] -= lambda[i][0];
			A[i][i] += lambda[i][0];
			A[i][i] -= lambda[i][1];
			A[i][i + 1] -= lambda[i][2];
			A[i][i] += lambda[i][2];
		}
		A[99][99] += 1.0 / dt;
		B[99] += X[99] / dt;
		//Tridiagonal matrix algorithm
		C[0] = A[i][i + 1] / A[i][i];
		for (i = 1; i < 99; i++)
			C[i] = A[i][i + 1] / (A[i][i] - A[i][i - 1] * C[i - 1]);
		D[0] = B[0] / A[i][i];
		for (i = 1; i < 100; i++)
			D[i] = (B[i] - A[i][i - 1] * D[i - 1]) / (A[i][i] - A[i][i - 1] * C[i - 1]);
		X[99] = D[99];
		for (i = 98; i >= 0; i--)
			X[i] = D[i] - C[i] * X[i + 1];
		
		sprintf(file_name, "tmp/trans_A_%d.dat", t);
		f = fopen(file_name, "w");
		for (i = 0; i < 100; i++)
			for (j = 0; j < 100; j++)
				fprintf(f, "%10.10lf\n", A[i][j]);
		fclose(f);
	}
	return 0;
}
