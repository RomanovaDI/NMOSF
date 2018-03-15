#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int main(void)
{
	double A[100][100], B[100], X[100];
	double dt = 0.1;
	double dx = 1;
	int i, t, j;
	double end_T = 10;
	FILE *f;
	char file_name[20 + i];
	memset(X, 0, 100 * sizeof(double));
	for (i = 30; i < 40; i++)
		X[i] = 1;
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
				A[i][i-1] += 1.0 / dx;
				A[i][i] += -1.0 / dx;
			}
		}
		for (i = 0; i < 100; i++) {
			if (i == 0)
				X[i] = B[i] / A[i][i];
			else
				X[i] = (B[i] - A[i][i-1] * X[i-1]) / A[i][i];
		}
		*/
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
		sprintf(file_name, "tmp/trans_A_%d.dat", t);
		f = fopen(file_name, "w");
		for (i = 0; i < 100; i++)
			for (j = 0; j < 100; j++)
				fprintf(f, "%10.10lf\n", A[i][j]);
		fclose(f);
	}
	return 0;
}
