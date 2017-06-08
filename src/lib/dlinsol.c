
/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */
#include "slu_ddefs.h"

int main(int argc, char *argv[])
{
    SuperMatrix A;
    NCformat *Astore;
    double   *a;
    int      *asub, *xa;
    int      *perm_c; /* column permutation vector */
    int      *perm_r; /* row permutations from partial pivoting */
    SuperMatrix L;      /* factor L */
    SCformat *Lstore;
    SuperMatrix U;      /* factor U */
    NCformat *Ustore;
    SuperMatrix B;
    int      nrhs, ldx, info, m, n, nnz;
    double   *xact, *rhs;
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
    //dreadhb(&m, &n, &nnz, &a, &asub, &xa);
#if 0
	m = n = 5;
	nnz = 11;
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
	a[0] = 1;
	a[1] = 1;
	a[2] = 1;
	a[3] = 2;
	a[4] = 1;
	a[5] = 1;
	a[6] = 3;
	a[7] = 1;
	a[8] = 4;
	a[9] = 1;
	a[10] = 5;
	asub[0] = 0;
	asub[1] = 2;
	asub[2] = 4;
	asub[3] = 1;
	asub[4] = 3;
	asub[5] = 1;
	asub[6] = 2;
	asub[7] = 2;
	asub[8] = 3;
	asub[9] = 0;
	asub[10] = 4;
	xa[0] = 0;
	xa[1] = 3;
	xa[2] = 5;
	xa[3] = 7;
	xa[4] = 9;
	xa[5] = 11;


    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");

	rhs[0] = 6;
	rhs[1] = 7;
	rhs[2] = 14;
	rhs[3] = 18;
	rhs[4] = 26;

    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);

	xact[0] = 1;
	xact[1] = 2;
	xact[2] = 3;
	xact[3] = 4;
	xact[4] = 5;
#endif

#if 1
	m = n = 5;
	nnz = 5;
    if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
	a[0] = 1;
	a[1] = 1;
	a[2] = 1;
	a[3] = 1;
	a[4] = 1;
	asub[0] = 0;
	asub[1] = 1;
	asub[2] = 2;
	asub[3] = 3;
	asub[4] = 4;
	xa[0] = 0;
	xa[1] = 1;
	xa[2] = 2;
	xa[3] = 3;
	xa[4] = 4;
	xa[5] = 5;


    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    Astore = A.Store;
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
    
    nrhs   = 1;
    if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");

	rhs[0] = 1;
	rhs[1] = 2;
	rhs[2] = 3;
	rhs[3] = 4;
	rhs[4] = 5;

    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);
    xact = doubleMalloc(n * nrhs);

	xact[0] = 1;
	xact[1] = 2;
	xact[2] = 3;
	xact[3] = 4;
	xact[4] = 5;
#endif

    ldx = n;
//    dGenXtrue(n, nrhs, xact, ldx);
//    dFillRHS(options.Trans, nrhs, xact, ldx, &A, &B);

    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    if ( info == 0 ) {

	/* This is how you could access the solution matrix. */
        double *sol = (double*) ((DNformat*) B.Store)->nzval; 

	 /* Compute the infinity norm of the error. */
	dinf_norm_error(nrhs, &B, xact);

	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
    	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
    	printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);
	
	dQuerySpace(&L, &U, &mem_usage);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	
    } else {
	printf("dgssv() error returns INFO= %d\n", info);
	if ( info <= n ) { /* factorization completes */
	    dQuerySpace(&L, &U, &mem_usage);
	    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
	}
    }

	dPrint_CompCol_Matrix("A", &A);
	dPrint_Dense_Matrix("B", &B);
	dPrint_CompCol_Matrix("U", &U);
	dPrint_SuperNode_Matrix("L", &L);
	print_int_vec("\nperm_r", m, perm_r);


    if ( options.PrintStat ) StatPrint(&stat);
    StatFree(&stat);

    SUPERLU_FREE (rhs);
    SUPERLU_FREE (xact);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC("Exit main()");
#endif
}

