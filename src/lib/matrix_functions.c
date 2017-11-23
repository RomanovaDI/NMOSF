#include "init_data.h"
#include "utils.h"
#include "matrix_functions.h"
#include "slu_ddefs.h"
#include <mpi.h>
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
	if ((f = fopen("tmp/A_csr.txt","w")) == NULL) {
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
	if ((f = fopen("tmp/A_B.txt","w")) == NULL) {
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
	if ((f = fopen("tmp/A_pattern.dat","w")) == NULL) {
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
		while (!((i >= I->Aiptr_csr[j]) && (i < I->Aiptr_csr[j + 1])))
			j++;
		//A_pattern[j * I->system_dimension + I->Ajptr_csr[i]] = 1;
		A_pattern[(I->system_dimension - 1 - j) * I->system_dimension + I->Ajptr_csr[i]] = 1;
	}
	if ((f = fopen("tmp/A_pattern_matrix.dat","w")) == NULL) {
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
	if ((f = fopen("tmp/A_elem_jprt.txt","w")) == NULL) {
		printf("error openning file");
		return;
	}
	fprintf(f, "#%10d\n", I->system_dimension);
	for (i = 0; i < I->non_zero_elem; i++) {
		fprintf(f, "%10d\t%20.10lf\t%10d\n", i, I->Aelem_csr[i], I->Ajptr_csr[i]);
	}
	fclose(f);
	if ((f = fopen("tmp/A_iptr.txt","w")) == NULL) {
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

int SuperLU_solver(in *I)
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
		printf("Matrix has no solution\n");
		return 1;
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

int csr_to_csc(double *csr_a, int *csr_row_ptr, int *csr_col_ind, int nonzeros, int dim, double *csc_a, int *csc_col_ptr, int *csc_row_ind)
{
	int i, j, a_ind = 0, row_ind;
	for (j = 0; j < dim; j++) {
		csc_col_ptr[j] = a_ind;
		for (i = 0; i < nonzeros; i++) {
			if (csr_col_ind[i] == j) {
				for (row_ind = 0; row_ind < dim; row_ind++)
					if ((i >= csr_row_ptr[row_ind]) && (i < csr_row_ptr[row_ind + 1]))
						break;
				csc_a[a_ind] = csr_a[i];
				csc_row_ind[a_ind++] = row_ind;
			}
		}
	}
	return 0;
}

int csr_multiply_by_csc_eq_csr(double *csr_a, int *csr_row_ptr, int *csr_col_ind, int csr_nonzeros, double *csc_a, int *csc_col_ptr, int *csc_row_ind, int csc_nonzeros, int dim)
{
	int i;
	return 0;
}

int handmade_solver(in *I)
{
	int i, j, k, *tmp_csc_row_ind, *tmp_csc_col_ptr;
	double *tmp_csc_a, *tmp_up_tr_matr, cosa, sina, diag, coef1, coef2;
	if ((tmp_csc_a = (double *) malloc(I->non_zero_elem * sizeof(double))) == NULL) {
		printf("Memory error in function %s\n", __func__);
		return 1;
	}
	if ((tmp_csc_row_ind = (int *) malloc(I->non_zero_elem * sizeof(int))) == NULL) {
		printf("Memory error in function %s\n", __func__);
		return 1;
	}
	if ((tmp_csc_col_ptr = (int *) malloc((I->system_dimension + 1) * sizeof(int))) == NULL) {
		printf("Memory error in function %s\n", __func__);
		return 1;
	}
	if ((tmp_up_tr_matr = (double *) malloc((I->system_dimension + 1) * I->system_dimension * 0.5 * sizeof(double))) == NULL) {
		printf("Memory error in function %s\n", __func__);
		return 1;
	}
	if (csr_to_csc(I->Aelem_csr, I->Aiptr_csr, I->Ajptr_csr, I->non_zero_elem, I->system_dimension, tmp_csc_a, tmp_csc_col_ptr, tmp_csc_row_ind)) return 1;
	for (i = 0; i < I->system_dimension; i++) {
		for (j = tmp_csc_col_ptr[i]; j < tmp_csc_col_ptr[i + 1]; j++) {
			diag = 0;
			if (tmp_csc_row_ind[j] == i) diag = tmp_csc_a[j];
			if (tmp_csc_row_ind[j] < i) {
				cosa = diag / (sqrt(diag * diag + tmp_csc_a[j] * tmp_csc_a[j]));
				sina = - tmp_csc_a[j] / (sqrt(diag * diag + tmp_csc_a[j] * tmp_csc_a[j]));
				for (k = tmp_csc_col_ptr[i]; k < tmp_csc_col_ptr[i + 1]; k++) {
					coef1 = coef2 = 0;
					if (tmp_csc_row_ind[k] == i)
						coef1 = tmp_csc_a[k];
					if (tmp_csc_row_ind[k] == j)
						coef2 = tmp_csc_a[k];
				}
			}
		}
	}
	return 0;
}

int solve_matrix(in *I)
{
//	if (SuperLU_solver(I)) return 1;
	if (handmade_solver(I)) return 1;
	return 0;
}

int solve_test_matrix(void)
{
	printf("Solving matrix\n");
	SuperMatrix A_csr;
	NCformat *Astore;
	int nonzeros = 3, systemdim = 3;
	double a_csr[nonzeros], b[systemdim];
	int a_iptr[systemdim + 1], a_jptr[nonzeros];
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

	for (i = 0; i < nonzeros; i++)
		a_csr[i] = 1;
	a_jptr[0] = 0;
	a_jptr[1] = 1;
	a_jptr[2] = 2;
//	a_jptr[3] = 3;
//	a_jptr[4] = 2;
//	a_jptr[5] = 0;
//	a_jptr[6] = 1;
//	a_jptr[7] = 2;
	a_iptr[0] = 0;
	a_iptr[1] = 1;
	a_iptr[2] = 2;
	a_iptr[3] = 3;
	b[0] = b[1] = 3;
	b[2] = 4;
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
	dCreate_CompCol_Matrix(&A_csr, systemdim, systemdim, nonzeros, a_csr, a_jptr, a_iptr, SLU_NR, SLU_D, SLU_GE);
	Astore = A_csr.Store;
	printf("Dimension %dx%d; # nonzeros %d\n", A_csr.nrow, A_csr.ncol, Astore->nnz);
	nrhs   = 1;
//	if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
	dCreate_Dense_Matrix(&B_csr, systemdim, nrhs, b, systemdim, SLU_DN, SLU_D, SLU_GE);
//	xact = doubleMalloc(n * nrhs);
//	ldx = n;
//	dGenXtrue(n, nrhs, xact, ldx);
//	dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);
	if ( !(perm_c = intMalloc(systemdim)) ) ABORT("Malloc fails for perm_c[].");
	if ( !(perm_r = intMalloc(systemdim)) ) ABORT("Malloc fails for perm_r[].");

	/* Initialize the statistics variables. */
	StatInit(&stat);
	dgssv(&options, &A_csr, perm_c, perm_r, &L, &U, &B_csr, &stat, &info);
	if ( info == 0 ) {
		/* This is how you could access the solution matrix. */
		double *sol = (double*) ((DNformat*) B_csr.Store)->nzval; 
		for (i = 0; i < systemdim; i++) {
			b[i] = sol[i];
			printf("%lf\n", sol[i]);
		}
		/* Compute the infinity norm of the error. */
//		dinf_norm_error(nrhs, &B, xact);
		dinf_norm_error(nrhs, &B_csr, sol);
		Lstore = (SCformat *) L.Store;
		Ustore = (NCformat *) U.Store;
		printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
		printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
		printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - systemdim);
		printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - systemdim)/nonzeros);
		dQuerySpace(&L, &U, &mem_usage);
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	} else {
		printf("dgssv() error returns INFO= %d\n", info);
		if ( info <= 3 ) { /* factorization completes */
			dQuerySpace(&L, &U, &mem_usage);
			printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
					mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		}
	}
	if ( options.PrintStat ) StatPrint(&stat);
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

