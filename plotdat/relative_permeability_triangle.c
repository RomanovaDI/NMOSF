#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void)
{
	int i;
	double s[3], sr[3];
	FILE *f[3];
	f[0] = fopen("plotdat/relative_permeability_triangle_water.dat", "w");
	f[1] = fopen("plotdat/relative_permeability_triangle_oil.dat", "w");
	f[2] = fopen("plotdat/relative_permeability_triangle_gas.dat", "w");
	sr[0] = 0.15;
	sr[1] = 0.2;
	sr[2] = 0;
	double k[3], kr[6], S;
	double k_max[3], k_min[3];
	i = 0;
	double step = 0.01;
	for (s[0] = sr[0] + step; s[0] < 1 - sr[1] - step - sr[2] - step; s[0] += step) {
		for (s[1] = sr[1] + step; s[1] < 1 - s[0] - sr[2] - step; s[1] += step) {
			for (s[2] = sr[2] + 0.01; s[2] < 1 - s[0] - s[1]; s[2] += 0.01) {
				S = (s[1] - sr[1]) / (1 - sr[1] - sr[0]);
				kr[0] = pow(S, 4);
				kr[1] = pow(1 - S, 2) * (1 - S * S);
				S = (s[0] - sr[0]) / (1 - sr[0] - sr[2]);
				kr[2] = pow(S, 4);
				kr[3] = pow(1 - S, 2) * (1 - S * S);
				S = (s[0] - sr[0]) / (1 - sr[0] - sr[1]);
				kr[4] = pow(S, 4);
				kr[5] = pow(1 - S, 2) * (1 - S * S);
				k[0] = ((s[1] - sr[1]) * kr[2] + (s[2] - sr[2]) * kr[4]) / (s[1] - sr[1] + s[2] - sr[2]);
				k[1] = ((s[0] - sr[0]) * kr[3] + (s[2] - sr[2]) * kr[0]) / (s[0] - sr[0] + s[2] - sr[2]);
				k[2] = ((s[0] - sr[0]) * kr[5] + (s[1] - sr[1]) * kr[1]) / (s[0] - sr[0] + s[1] - sr[1]);
				if (i == 0)
					for (i = 0; i < 3; i++)
						k_max[i] = k_min[i] = k[i];
				for (i = 0; i < 3; i++) {
					if (k_max[i] < k[i])
						k_max[i] = k[i];
					if (k_min[i] > k[i])
						k_min[i] = k[i];
					fprintf(f[i], "%10.10lf\t%10.10lf\t%10.10lf\t%10.10lf\n", s[0], s[1], s[2], k[i]);
				}
			}
		}
	}
	for (i = 0; i < 3; i++) {
		printf("Phase %d: k_min = %10.10lf, k_max = %10.10lf\n", i, k_min[i], k_max[i]);
		fclose(f[i]);
	}
	return 0;
}
