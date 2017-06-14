#include "init_data.h"
#include "utils.h"
#include "matrix_functions.h"
#include "slu_ddefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

int write_B_to_B_prev(in *I)
{
	int i, j, k, p;
	for (p = 0; p < I->num_parameters; p++) {
		for (k = 0; k < I->nz; k++) {
			for (i = 0; i < I->nx; i++) {
				for (j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						I->B_prev[B_IND(I, p, i, j, k)] = I->B[A_IND(I, p, i, j, k)];
				}
			}
		}
	}
	return 0;
}

int write_B_prev_to_B(in *I)
{
	int i, j, k, p;
	for (p = 0; p < I->num_parameters; p++) {
		for (k = 0; k < I->nz; k++) {
			for (i = 0; i < I->nx; i++) {
				for (j = 0; j < I->ny; j++) {
					if (I->ind_cell_multipl[i * I->ny + j] != -1)
						I->B[A_IND(I, p, i, j, k)] = I->B_prev[B_IND(I, p, i, j, k)];
				}
			}
		}
	}
	return 0;
}

void print_A_csr(in *I)
{
	printf("Print matrix A in CSR format function\n");
	FILE *f;
	int i, j, k, fl_tmp;
	if ((f = fopen("A_csr.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", I->system_dimension);
	double max = abs(I->Aelem_csr[0]);
	for (i = 0; i < I->non_zero_elem; i++) {
		if (max < abs(I->Aelem_csr[i]))
			max = abs(I->Aelem_csr[i]);
		j = 0;
		while (!((i >= I->Aiptr_csr[j]) && (i < I->Aiptr_csr[j + 1])))
			j++;
		fprintf(f, "%20.10lf\t%10d\t%10d\n", I->Aelem_csr[i], j, I->Ajptr_csr[i]);
	}
	fclose(f);
	if ((f = fopen("A_B.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	for (i = 0; i < I->system_dimension; i++) {
		for (j = 0; j < I->system_dimension; j++) {
			fl_tmp = 1;
			for (k = I->Aiptr_csr[i]; k < I->Aiptr_csr[i + 1]; k++) {
				if (I->Ajptr_csr[k] == j) {
					fprintf(f, "%20.10lf\t", I->Aelem_csr[k]);
					fl_tmp = 0;
					break;
				}
				if (fl_tmp) {
					fprintf(f, "%20.10lf\t", 0.0);
				}
			}
		}
		fprintf(f, "\t\t%20.10lf\n", I->B[i]);
	}
	fclose(f);
	if ((f = fopen("A_pattern.dat","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", I->system_dimension);
	for (i = 0; i < I->non_zero_elem; i++) {
		j = 0;
		while (!((i >= I->Aiptr_csr[j]) && (i <= I->Aiptr_csr[j + 1])))
			j++;
		//fprintf(f, "%20.10lf\t%10d\t%10d\n", I->Aelem_csr[i] / max, I->system_dimension - 1 - j, I->Ajptr_csr[i]);
		fprintf(f, "%20.10lf\t%10d\t%10d\n", 1., I->system_dimension - 1 - j, I->Ajptr_csr[i]);
	}
	fclose(f);
	int *A_pattern = (int *) malloc(I->system_dimension * I->system_dimension * sizeof(int));
	memset(A_pattern, 0, I->system_dimension * I->system_dimension * sizeof(int));
	for (i = 0; i < I->non_zero_elem; i++) {
		j = 0;
		while (!((i >= I->Aiptr_csr[j]) && (i <= I->Aiptr_csr[j + 1])))
			j++;
		A_pattern[(I->system_dimension - 1 - j) * I->system_dimension + I->Ajptr_csr[i]] = 1;
	}
	if ((f = fopen("A_pattern_matrix.dat","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", I->system_dimension);
	for (i = 0; i < I->system_dimension; i++) {
		for (j = 0; j < I->system_dimension; j++) {
			fprintf(f, "%d\t", A_pattern[i * I->system_dimension + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	if ((f = fopen("A_elem_jprt.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", I->system_dimension);
	for (i = 0; i < I->non_zero_elem; i++) {
		fprintf(f, "%10d\t%20.10lf\t%10d\n", i, I->Aelem_csr[i], I->Ajptr_csr[i]);
	}
	fclose(f);
	if ((f = fopen("A_iptr.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", I->system_dimension);
	for (i = 0; i < I->system_dimension; i++) {
		fprintf(f, "%10d\n", I->Aiptr_csr[i]);
	}
	fclose(f);
	return;
}

int solve_matrix(in *I)
{
	printf("Solving matrix\n");
	SuperMatrix A_csr;
	NCformat *Astore;
//	double   *a;
//	int      *asub, *xa;
	int      *perm_c; /* column permutation vector */
	int      *perm_r; /* row permutations from partial pivoting */
	SuperMatrix L;      /* factor L */
	SCformat *Lstore;
	SuperMatrix U;      /* factor U */
	NCformat *Ustore;
	SuperMatrix B_csr;
//	int      nrhs, ldx, info, m, n, nnz;
	int      info, i, j, nrhs;
//	double   *xact, *rhs;
//	double   *xact;
	mem_usage_t   mem_usage;
	superlu_options_t options;
	SuperLUStat_t stat;

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Enter main()");
#endif

	/* Set the default input options:
	options.Fact = DOFACT;
	options.Equil = YES;
	options.ColPerm = COLAMD;
	options.DiagPivotThresh = 1.0;
	options.Trans = NOTRANS;
	options.IterRefine = NOREFINE;
	options.SymmetricMode = NO;
	options.PivotGrowth = NO;
	options.ConditionNumber = NO;
	options.PrintStat = YES;
	 */

	set_default_options(&options);

	/* Read the matrix in Harwell-Boeing format. */
//	dreadhb(&m, &n, &nnz, &a, &asub, &xa);
	dCreate_CompCol_Matrix(&A_csr, I->system_dimension, I->system_dimension, I->non_zero_elem, I->Aelem_csr, I->Ajptr_csr, I->Aiptr_csr, SLU_NR, SLU_D, SLU_GE);
	Astore = A_csr.Store;
	printf("Dimension %dx%d; # nonzeros %d\n", A_csr.nrow, A_csr.ncol, Astore->nnz);
	nrhs   = 1;
//	if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
	dCreate_Dense_Matrix(&B_csr, I->system_dimension, nrhs, I->B, I->system_dimension, SLU_DN, SLU_D, SLU_GE);
//	xact = doubleMalloc(n * nrhs);
//	ldx = n;
//	dGenXtrue(n, nrhs, xact, ldx);
//	dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
	if ( !(perm_c = intMalloc(I->system_dimension)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(perm_r = intMalloc(I->system_dimension)) ) ABORT("Malloc fails for perm_r[].");

	/* Initialize the statistics variables. */
	StatInit(&stat);
	dgssv(&options, &A_csr, perm_c, perm_r, &L, &U, &B_csr, &stat, &info);
	if ( info == 0 ) {
		/* This is how you could access the solution matrix. */
		double *sol = (double*) ((DNformat*) B_csr.Store)->nzval; 
		for (i = 0; i < I->n_cells_multipl * I->nz * I->num_parameters; i++) {
			I->B[i] = sol[i];
		}
		/* Compute the infinity norm of the error. */
//		dinf_norm_error(nrhs, &B, xact);
		dinf_norm_error(nrhs, &B_csr, sol);
		Lstore = (SCformat *) L.Store;
		Ustore = (NCformat *) U.Store;
		printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
		printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
		printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - I->system_dimension);
		printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - I->system_dimension)/I->non_zero_elem);
		dQuerySpace(&L, &U, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	} else {
		printf("dgssv() error returns INFO= %d\n", info);
		if ( info <= I->system_dimension ) { /* factorization completes */
			dQuerySpace(&L, &U, &mem_usage);
			printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
					mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		}
	}
	if ( options.PrintStat ) StatPrint(&stat);
//	for (i = 0; i < n_cells_multipl * nz; i++) {
//		velocity[i * 3 + 0] = sol[i * num_parameters + 0];
//		velocity[i * 3 + 1] = sol[i * num_parameters + 1];
//		velocity[i * 3 + 3] = sol[i * num_parameters + 2];
//		phase_fraction[i] = sol[i * num_parameters + 3];
//		pressure[i] = sol[i * num_parameters + 4];
//	}
	printf("Matrix solved\n");
	
	StatFree(&stat);
//	SUPERLU_FREE (rhs);
//	SUPERLU_FREE (xact);
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);
//	Destroy_CompCol_Matrix(&A_csr);
//	Destroy_SuperMatrix_Store(&B_csr);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC("Exit main()");
#endif
	printf("End of function solve matrix\n");
	return 0;
}
