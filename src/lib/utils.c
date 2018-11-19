#include "init_data.h"
#include "utils.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

int A_IND(in *I, int p, int i, int j, int k)
{
	return (I->num_parameters * (k * I->n_cells_multipl + I->ind_cell_multipl[i * I->ny + j]) + p); // index in vector of parameters without boundaries
}

int GL_A_IND(in *I, int p, int i, int j, int k)
{
#if DEBUG==4
	printf("I->gl_n_cells_multipl = %d\n", I->gl_n_cells_multipl);
	printf("Processor %d PID %d: i = %d, j = %d, k = %d, p = %d, I->gl_n_cells_multipl = %d\n", I->my_rank, getpid(), i, j, k, p, I->gl_n_cells_multipl);
	printf("Processor %d PID %d: i = %d, j = %d, k = %d, p = %d, I->gl_ind_cell_multipl[i * I->gl_ny + j] = %d\n", I->my_rank, getpid(), i, j, k, p, I->gl_ind_cell_multipl[i * I->gl_ny + j]);
	printf("Processor %d PID %d: i = %d, j = %d, k = %d, p = %d, I->gl_n_cells_multipl = %d, I->gl_ind_cell_multipl[i * I->gl_ny + j] = %d\n", I->my_rank, getpid(), i, j, k, p, I->gl_n_cells_multipl, I->gl_ind_cell_multipl[i * I->gl_ny + j]);
#endif
	return (I->num_parameters * (k * I->gl_n_cells_multipl + I->gl_ind_cell_multipl[i * I->gl_ny + j]) + p);
}

int B_IND(in *I, int p, int i, int j, int k) // index in vector of parameters with boundaries
{
	return (I->num_parameters * ((k + I->stencil_size) * I->n_boundary_cells + I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size]) + p);
}

int PARAM_IND(in *I, int p, int i, int j, int k)
{
	return (I->dependent_variables * ((k + I->stencil_size) * I->n_boundary_cells + I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size]) + p);
}

int POR_IND(in *I, int i, int j, int k) // porousness index
{
	return (k + I->stencil_size) * I->n_boundary_cells + I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size];
}

int INT_NEI_IND(in *I, int p, int i, int j, int k) // interior neighbor cell index
{
	if (internal_cell(I, i, j, k))
		return B_IND(I, p, i, j, k);
	else if (internal_cell(I, i + 1, j, k))
		return B_IND(I, p, i + 1, j, k);
	else if (internal_cell(I, i - 1, j, k))
		return B_IND(I, p, i - 1, j, k);
	else if (internal_cell(I, i, j + 1, k))
		return B_IND(I, p, i, j + 1, k);
	else if (internal_cell(I, i, j - 1, k))
		return B_IND(I, p, i, j - 1, k);
	else if (internal_cell(I, i, j, k + 1))
		return B_IND(I, p, i, j, k + 1);
	else if (internal_cell(I, i, j, k - 1))
		return B_IND(I, p, i, j, k - 1);
	else if (internal_cell(I, i + 1, j + 1, k))
		return B_IND(I, p, i + 1, j + 1, k);
	else if (internal_cell(I, i - 1, j + 1, k))
		return B_IND(I, p, i - 1, j + 1, k);
	else if (internal_cell(I, i + 1, j - 1, k))
		return B_IND(I, p, i + 1, j - 1, k);
	else if (internal_cell(I, i - 1, j - 1, k))
		return B_IND(I, p, i - 1, j - 1, k);
	else if (internal_cell(I, i + 1, j, k + 1))
		return B_IND(I, p, i + 1, j, k + 1);
	else if (internal_cell(I, i - 1, j, k + 1))
		return B_IND(I, p, i - 1, j, k + 1);
	else if (internal_cell(I, i + 1, j, k - 1))
		return B_IND(I, p, i + 1, j, k - 1);
	else if (internal_cell(I, i - 1, j, k - 1))
		return B_IND(I, p, i - 1, j, k - 1);
	else if (internal_cell(I, i, j + 1, k + 1))
		return B_IND(I, p, i, j + 1, k + 1);
	else if (internal_cell(I, i, j - 1, k + 1))
		return B_IND(I, p, i, j - 1, k + 1);
	else if (internal_cell(I, i, j + 1, k - 1))
		return B_IND(I, p, i, j + 1, k - 1);
	else if (internal_cell(I, i, j - 1, k - 1))
		return B_IND(I, p, i, j - 1, k - 1);
	else if (internal_cell(I, i + 1, j + 1, k + 1))
		return B_IND(I, p, i + 1, j + 1, k + 1);
	else if (internal_cell(I, i - 1, j + 1, k + 1))
		return B_IND(I, p, i - 1, j + 1, k + 1);
	else if (internal_cell(I, i + 1, j - 1, k + 1))
		return B_IND(I, p, i + 1, j - 1, k + 1);
	else if (internal_cell(I, i + 1, j + 1, k - 1))
		return B_IND(I, p, i + 1, j + 1, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k + 1))
		return B_IND(I, p, i - 1, j - 1, k + 1);
	else if (internal_cell(I, i - 1, j + 1, k - 1))
		return B_IND(I, p, i - 1, j + 1, k - 1);
	else if (internal_cell(I, i + 1, j - 1, k - 1))
		return B_IND(I, p, i + 1, j - 1, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k - 1))
		return B_IND(I, p, i - 1, j - 1, k - 1);
	else if (internal_cell(I, i + 2, j, k))
		return B_IND(I, p, i + 2, j, k);
	else if (internal_cell(I, i - 2, j, k))
		return B_IND(I, p, i - 2, j, k);
	else if (internal_cell(I, i, j + 2, k))
		return B_IND(I, p, i, j + 2, k);
	else if (internal_cell(I, i, j - 2, k))
		return B_IND(I, p, i, j - 2, k);
	else if (internal_cell(I, i, j, k + 2))
		return B_IND(I, p, i, j, k + 2);
	else if (internal_cell(I, i, j, k - 2))
		return B_IND(I, p, i, j, k - 2);
	else if (internal_cell(I, i + 2, j + 1, k))
		return B_IND(I, p, i + 2, j + 1, k);
	else if (internal_cell(I, i + 2, j, k + 1))
		return B_IND(I, p, i + 2, j, k + 1);
	else if (internal_cell(I, i + 1, j + 2, k))
		return B_IND(I, p, i + 1, j + 2, k);
	else if (internal_cell(I, i, j + 2, k + 1))
		return B_IND(I, p, i, j + 2, k + 1);
	else if (internal_cell(I, i + 1, j, k + 2))
		return B_IND(I, p, i + 1, j, k + 2);
	else if (internal_cell(I, i, j + 1, k + 2))
		return B_IND(I, p, i, j + 1, k + 2);
	else if (internal_cell(I, i - 2, j + 1, k))
		return B_IND(I, p, i - 2, j + 1, k);
	else if (internal_cell(I, i - 2, j, k + 1))
		return B_IND(I, p, i - 2, j, k + 1);
	else if (internal_cell(I, i + 1, j - 2, k))
		return B_IND(I, p, i + 1, j - 2, k);
	else if (internal_cell(I, i, j - 2, k + 1))
		return B_IND(I, p, i, j - 2, k + 1);
	else if (internal_cell(I, i + 1, j, k - 2))
		return B_IND(I, p, i + 1, j, k - 2);
	else if (internal_cell(I, i, j + 1, k - 2))
		return B_IND(I, p, i, j + 1, k - 2);
	else if (internal_cell(I, i + 2, j - 1, k))
		return B_IND(I, p, i + 2, j - 1, k);
	else if (internal_cell(I, i + 2, j, k - 1))
		return B_IND(I, p, i + 2, j, k - 1);
	else if (internal_cell(I, i - 1, j + 2, k))
		return B_IND(I, p, i - 1, j + 2, k);
	else if (internal_cell(I, i, j + 2, k - 1))
		return B_IND(I, p, i, j + 2, k - 1);
	else if (internal_cell(I, i - 1, j, k + 2))
		return B_IND(I, p, i - 1, j, k + 2);
	else if (internal_cell(I, i, j - 1, k + 2))
		return B_IND(I, p, i, j - 1, k + 2);
	else if (internal_cell(I, i - 2, j - 1, k))
		return B_IND(I, p, i - 2, j - 1, k);
	else if (internal_cell(I, i - 2, j, k - 1))
		return B_IND(I, p, i - 2, j, k - 1);
	else if (internal_cell(I, i - 1, j - 2, k))
		return B_IND(I, p, i - 1, j - 2, k);
	else if (internal_cell(I, i, j - 2, k - 1))
		return B_IND(I, p, i, j - 2, k - 1);
	else if (internal_cell(I, i - 1, j, k - 2))
		return B_IND(I, p, i - 1, j, k - 2);
	else if (internal_cell(I, i, j - 1, k - 2))
		return B_IND(I, p, i, j - 1, k - 2);
	else if (internal_cell(I, i + 2, j + 2, k))
		return B_IND(I, p, i + 2, j + 2, k);
	else if (internal_cell(I, i - 2, j + 2, k))
		return B_IND(I, p, i - 2, j + 2, k);
	else if (internal_cell(I, i + 2, j - 2, k))
		return B_IND(I, p, i + 2, j - 2, k);
	else if (internal_cell(I, i - 2, j - 2, k))
		return B_IND(I, p, i - 2, j - 2, k);
	else if (internal_cell(I, i + 2, j, k + 2))
		return B_IND(I, p, i + 2, j, k + 2);
	else if (internal_cell(I, i - 2, j, k + 2))
		return B_IND(I, p, i - 2, j, k + 2);
	else if (internal_cell(I, i + 2, j, k - 2))
		return B_IND(I, p, i + 2, j, k - 2);
	else if (internal_cell(I, i - 2, j, k - 2))
		return B_IND(I, p, i - 2, j, k - 2);
	else if (internal_cell(I, i, j + 2, k + 2))
		return B_IND(I, p, i, j + 2, k + 2);
	else if (internal_cell(I, i, j - 2, k + 2))
		return B_IND(I, p, i, j - 2, k + 2);
	else if (internal_cell(I, i, j + 2, k - 2))
		return B_IND(I, p, i, j + 2, k - 2);
	else if (internal_cell(I, i, j - 2, k - 2))
		return B_IND(I, p, i, j - 2, k - 2);
	else if (internal_cell(I, i + 1, j + 1, k + 2))
		return B_IND(I, p, i + 1, j + 1, k + 2);
	else if (internal_cell(I, i + 1, j + 2, k + 1))
		return B_IND(I, p, i + 1, j + 2, k + 1);
	else if (internal_cell(I, i + 2, j + 1, k + 1))
		return B_IND(I, p, i + 2, j + 1, k + 1);
	else if (internal_cell(I, i - 2, j + 1, k + 1))
		return B_IND(I, p, i - 2, j + 1, k + 1);
	else if (internal_cell(I, i + 1, j - 2, k + 1))
		return B_IND(I, p, i + 1, j - 2, k + 1);
	else if (internal_cell(I, i + 1, j + 1, k - 2))
		return B_IND(I, p, i + 1, j + 1, k - 2);
	else if (internal_cell(I, i + 2, j - 1, k + 1))
		return B_IND(I, p, i + 2, j - 1, k + 1);
	else if (internal_cell(I, i + 2, j + 1, k - 1))
		return B_IND(I, p, i + 2, j + 1, k - 1);
	else if (internal_cell(I, i - 1, j + 2, k + 1))
		return B_IND(I, p, i - 1, j + 2, k + 1);
	else if (internal_cell(I, i + 1, j + 2, k - 1))
		return B_IND(I, p, i + 1, j + 2, k - 1);
	else if (internal_cell(I, i - 1, j + 1, k + 2))
		return B_IND(I, p, i - 1, j + 1, k + 2);
	else if (internal_cell(I, i + 1, j - 1, k + 2))
		return B_IND(I, p, i + 1, j - 1, k + 2);
	else if (internal_cell(I, i - 2, j - 1, k + 1))
		return B_IND(I, p, i - 2, j - 1, k + 1);
	else if (internal_cell(I, i - 2, j + 1, k - 1))
		return B_IND(I, p, i - 2, j + 1, k - 1);
	else if (internal_cell(I, i - 1, j - 2, k + 1))
		return B_IND(I, p, i - 1, j - 2, k + 1);
	else if (internal_cell(I, i + 1, j - 2, k - 1))
		return B_IND(I, p, i + 1, j - 2, k - 1);
	else if (internal_cell(I, i - 1, j + 1, k - 2))
		return B_IND(I, p, i - 1, j + 1, k - 2);
	else if (internal_cell(I, i + 1, j - 1, k - 2))
		return B_IND(I, p, i + 1, j - 1, k - 2);
	else if (internal_cell(I, i + 2, j - 1, k - 1))
		return B_IND(I, p, i + 2, j - 1, k - 1);
	else if (internal_cell(I, i - 1, j + 2, k - 1))
		return B_IND(I, p, i - 1, j + 2, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k + 2))
		return B_IND(I, p, i - 1, j - 1, k + 2);
	else if (internal_cell(I, i - 2, j - 1, k - 1))
		return B_IND(I, p, i - 2, j - 1, k - 1);
	else if (internal_cell(I, i - 1, j - 2, k - 1))
		return B_IND(I, p, i - 1, j - 2, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k - 2))
		return B_IND(I, p, i - 1, j - 1, k - 2);
	else if (internal_cell(I, i + 1, j + 2, k + 2))
		return B_IND(I, p, i + 1, j + 2, k + 2);
	else if (internal_cell(I, i + 2, j + 1, k + 2))
		return B_IND(I, p, i + 2, j + 1, k + 2);
	else if (internal_cell(I, i + 2, j + 2, k + 1))
		return B_IND(I, p, i + 2, j + 2, k + 1);
	else if (internal_cell(I, i - 1, j + 2, k + 2))
		return B_IND(I, p, i - 1, j + 2, k + 2);
	else if (internal_cell(I, i + 2, j - 1, k + 2))
		return B_IND(I, p, i + 2, j - 1, k + 2);
	else if (internal_cell(I, i + 2, j + 2, k - 1))
		return B_IND(I, p, i + 2, j + 2, k - 1);
	else if (internal_cell(I, i + 1, j - 2, k + 2))
		return B_IND(I, p, i + 1, j - 2, k + 2);
	else if (internal_cell(I, i + 1, j + 2, k - 2))
		return B_IND(I, p, i + 1, j + 2, k - 2);
	else if (internal_cell(I, i - 2, j + 1, k + 2))
		return B_IND(I, p, i - 2, j + 1, k + 2);
	else if (internal_cell(I, i + 2, j + 1, k - 2))
		return B_IND(I, p, i + 2, j + 1, k - 2);
	else if (internal_cell(I, i - 2, j + 2, k + 1))
		return B_IND(I, p, i - 2, j + 2, k + 1);
	else if (internal_cell(I, i + 2, j - 2, k + 1))
		return B_IND(I, p, i + 2, j - 2, k + 1);
	else if (internal_cell(I, i - 1, j - 2, k + 2))
		return B_IND(I, p, i - 1, j - 2, k + 2);
	else if (internal_cell(I, i - 1, j + 2, k - 2))
		return B_IND(I, p, i - 1, j + 2, k - 2);
	else if (internal_cell(I, i - 2, j - 1, k + 2))
		return B_IND(I, p, i - 2, j - 1, k + 2);
	else if (internal_cell(I, i + 2, j - 1, k - 2))
		return B_IND(I, p, i + 2, j - 1, k - 2);
	else if (internal_cell(I, i - 2, j + 2, k - 1))
		return B_IND(I, p, i - 2, j + 2, k - 1);
	else if (internal_cell(I, i + 2, j - 2, k - 1))
		return B_IND(I, p, i + 2, j - 2, k - 1);
	else if (internal_cell(I, i + 1, j - 2, k - 2))
		return B_IND(I, p, i + 1, j - 2, k - 2);
	else if (internal_cell(I, i - 2, j + 1, k - 2))
		return B_IND(I, p, i - 2, j + 1, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k + 1))
		return B_IND(I, p, i - 2, j - 2, k + 1);
	else if (internal_cell(I, i - 1, j - 2, k - 2))
		return B_IND(I, p, i - 1, j - 2, k - 2);
	else if (internal_cell(I, i - 2, j - 1, k - 2))
		return B_IND(I, p, i - 2, j - 1, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k - 1))
		return B_IND(I, p, i - 2, j - 2, k - 1);
	else if (internal_cell(I, i + 2, j + 2, k + 2))
		return B_IND(I, p, i + 2, j + 2, k + 2);
	else if (internal_cell(I, i - 2, j + 2, k + 2))
		return B_IND(I, p, i - 2, j + 2, k + 2);
	else if (internal_cell(I, i + 2, j - 2, k + 2))
		return B_IND(I, p, i + 2, j - 2, k + 2);
	else if (internal_cell(I, i + 2, j + 2, k - 2))
		return B_IND(I, p, i + 2, j + 2, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k + 2))
		return B_IND(I, p, i - 2, j - 2, k + 2);
	else if (internal_cell(I, i - 2, j + 2, k - 2))
		return B_IND(I, p, i - 2, j + 2, k - 2);
	else if (internal_cell(I, i + 2, j - 2, k - 2))
		return B_IND(I, p, i + 2, j - 2, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k - 2))
		return B_IND(I, p, i - 2, j - 2, k - 2);
}

int POR_INT_NEI_IND(in *I, int i, int j, int k) // porousness interior neighbor cell index
{
	if (internal_cell(I, i, j, k))
		return POR_IND(I, i, j, k);
	else if (internal_cell(I, i + 1, j, k))
		return POR_IND(I, i + 1, j, k);
	else if (internal_cell(I, i - 1, j, k))
		return POR_IND(I, i - 1, j, k);
	else if (internal_cell(I, i, j + 1, k))
		return POR_IND(I, i, j + 1, k);
	else if (internal_cell(I, i, j - 1, k))
		return POR_IND(I, i, j - 1, k);
	else if (internal_cell(I, i, j, k + 1))
		return POR_IND(I, i, j, k + 1);
	else if (internal_cell(I, i, j, k - 1))
		return POR_IND(I, i, j, k - 1);
	else if (internal_cell(I, i + 1, j + 1, k))
		return POR_IND(I, i + 1, j + 1, k);
	else if (internal_cell(I, i - 1, j + 1, k))
		return POR_IND(I, i - 1, j + 1, k);
	else if (internal_cell(I, i + 1, j - 1, k))
		return POR_IND(I, i + 1, j - 1, k);
	else if (internal_cell(I, i - 1, j - 1, k))
		return POR_IND(I, i - 1, j - 1, k);
	else if (internal_cell(I, i + 1, j, k + 1))
		return POR_IND(I, i + 1, j, k + 1);
	else if (internal_cell(I, i - 1, j, k + 1))
		return POR_IND(I, i - 1, j, k + 1);
	else if (internal_cell(I, i + 1, j, k - 1))
		return POR_IND(I, i + 1, j, k - 1);
	else if (internal_cell(I, i - 1, j, k - 1))
		return POR_IND(I, i - 1, j, k - 1);
	else if (internal_cell(I, i, j + 1, k + 1))
		return POR_IND(I, i, j + 1, k + 1);
	else if (internal_cell(I, i, j - 1, k + 1))
		return POR_IND(I, i, j - 1, k + 1);
	else if (internal_cell(I, i, j + 1, k - 1))
		return POR_IND(I, i, j + 1, k - 1);
	else if (internal_cell(I, i, j - 1, k - 1))
		return POR_IND(I, i, j - 1, k - 1);
	else if (internal_cell(I, i + 1, j + 1, k + 1))
		return POR_IND(I, i + 1, j + 1, k + 1);
	else if (internal_cell(I, i - 1, j + 1, k + 1))
		return POR_IND(I, i - 1, j + 1, k + 1);
	else if (internal_cell(I, i + 1, j - 1, k + 1))
		return POR_IND(I, i + 1, j - 1, k + 1);
	else if (internal_cell(I, i + 1, j + 1, k - 1))
		return POR_IND(I, i + 1, j + 1, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k + 1))
		return POR_IND(I, i - 1, j - 1, k + 1);
	else if (internal_cell(I, i - 1, j + 1, k - 1))
		return POR_IND(I, i - 1, j + 1, k - 1);
	else if (internal_cell(I, i + 1, j - 1, k - 1))
		return POR_IND(I, i + 1, j - 1, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k - 1))
		return POR_IND(I, i - 1, j - 1, k - 1);
	else if (internal_cell(I, i + 2, j, k))
		return POR_IND(I, i + 2, j, k);
	else if (internal_cell(I, i - 2, j, k))
		return POR_IND(I, i - 2, j, k);
	else if (internal_cell(I, i, j + 2, k))
		return POR_IND(I, i, j + 2, k);
	else if (internal_cell(I, i, j - 2, k))
		return POR_IND(I, i, j - 2, k);
	else if (internal_cell(I, i, j, k + 2))
		return POR_IND(I, i, j, k + 2);
	else if (internal_cell(I, i, j, k - 2))
		return POR_IND(I, i, j, k - 2);
	else if (internal_cell(I, i + 2, j + 1, k))
		return POR_IND(I, i + 2, j + 1, k);
	else if (internal_cell(I, i + 2, j, k + 1))
		return POR_IND(I, i + 2, j, k + 1);
	else if (internal_cell(I, i + 1, j + 2, k))
		return POR_IND(I, i + 1, j + 2, k);
	else if (internal_cell(I, i, j + 2, k + 1))
		return POR_IND(I, i, j + 2, k + 1);
	else if (internal_cell(I, i + 1, j, k + 2))
		return POR_IND(I, i + 1, j, k + 2);
	else if (internal_cell(I, i, j + 1, k + 2))
		return POR_IND(I, i, j + 1, k + 2);
	else if (internal_cell(I, i - 2, j + 1, k))
		return POR_IND(I, i - 2, j + 1, k);
	else if (internal_cell(I, i - 2, j, k + 1))
		return POR_IND(I, i - 2, j, k + 1);
	else if (internal_cell(I, i + 1, j - 2, k))
		return POR_IND(I, i + 1, j - 2, k);
	else if (internal_cell(I, i, j - 2, k + 1))
		return POR_IND(I, i, j - 2, k + 1);
	else if (internal_cell(I, i + 1, j, k - 2))
		return POR_IND(I, i + 1, j, k - 2);
	else if (internal_cell(I, i, j + 1, k - 2))
		return POR_IND(I, i, j + 1, k - 2);
	else if (internal_cell(I, i + 2, j - 1, k))
		return POR_IND(I, i + 2, j - 1, k);
	else if (internal_cell(I, i + 2, j, k - 1))
		return POR_IND(I, i + 2, j, k - 1);
	else if (internal_cell(I, i - 1, j + 2, k))
		return POR_IND(I, i - 1, j + 2, k);
	else if (internal_cell(I, i, j + 2, k - 1))
		return POR_IND(I, i, j + 2, k - 1);
	else if (internal_cell(I, i - 1, j, k + 2))
		return POR_IND(I, i - 1, j, k + 2);
	else if (internal_cell(I, i, j - 1, k + 2))
		return POR_IND(I, i, j - 1, k + 2);
	else if (internal_cell(I, i - 2, j - 1, k))
		return POR_IND(I, i - 2, j - 1, k);
	else if (internal_cell(I, i - 2, j, k - 1))
		return POR_IND(I, i - 2, j, k - 1);
	else if (internal_cell(I, i - 1, j - 2, k))
		return POR_IND(I, i - 1, j - 2, k);
	else if (internal_cell(I, i, j - 2, k - 1))
		return POR_IND(I, i, j - 2, k - 1);
	else if (internal_cell(I, i - 1, j, k - 2))
		return POR_IND(I, i - 1, j, k - 2);
	else if (internal_cell(I, i, j - 1, k - 2))
		return POR_IND(I, i, j - 1, k - 2);
	else if (internal_cell(I, i + 2, j + 2, k))
		return POR_IND(I, i + 2, j + 2, k);
	else if (internal_cell(I, i - 2, j + 2, k))
		return POR_IND(I, i - 2, j + 2, k);
	else if (internal_cell(I, i + 2, j - 2, k))
		return POR_IND(I, i + 2, j - 2, k);
	else if (internal_cell(I, i - 2, j - 2, k))
		return POR_IND(I, i - 2, j - 2, k);
	else if (internal_cell(I, i + 2, j, k + 2))
		return POR_IND(I, i + 2, j, k + 2);
	else if (internal_cell(I, i - 2, j, k + 2))
		return POR_IND(I, i - 2, j, k + 2);
	else if (internal_cell(I, i + 2, j, k - 2))
		return POR_IND(I, i + 2, j, k - 2);
	else if (internal_cell(I, i - 2, j, k - 2))
		return POR_IND(I, i - 2, j, k - 2);
	else if (internal_cell(I, i, j + 2, k + 2))
		return POR_IND(I, i, j + 2, k + 2);
	else if (internal_cell(I, i, j - 2, k + 2))
		return POR_IND(I, i, j - 2, k + 2);
	else if (internal_cell(I, i, j + 2, k - 2))
		return POR_IND(I, i, j + 2, k - 2);
	else if (internal_cell(I, i, j - 2, k - 2))
		return POR_IND(I, i, j - 2, k - 2);
	else if (internal_cell(I, i + 1, j + 1, k + 2))
		return POR_IND(I, i + 1, j + 1, k + 2);
	else if (internal_cell(I, i + 1, j + 2, k + 1))
		return POR_IND(I, i + 1, j + 2, k + 1);
	else if (internal_cell(I, i + 2, j + 1, k + 1))
		return POR_IND(I, i + 2, j + 1, k + 1);
	else if (internal_cell(I, i - 2, j + 1, k + 1))
		return POR_IND(I, i - 2, j + 1, k + 1);
	else if (internal_cell(I, i + 1, j - 2, k + 1))
		return POR_IND(I, i + 1, j - 2, k + 1);
	else if (internal_cell(I, i + 1, j + 1, k - 2))
		return POR_IND(I, i + 1, j + 1, k - 2);
	else if (internal_cell(I, i + 2, j - 1, k + 1))
		return POR_IND(I, i + 2, j - 1, k + 1);
	else if (internal_cell(I, i + 2, j + 1, k - 1))
		return POR_IND(I, i + 2, j + 1, k - 1);
	else if (internal_cell(I, i - 1, j + 2, k + 1))
		return POR_IND(I, i - 1, j + 2, k + 1);
	else if (internal_cell(I, i + 1, j + 2, k - 1))
		return POR_IND(I, i + 1, j + 2, k - 1);
	else if (internal_cell(I, i - 1, j + 1, k + 2))
		return POR_IND(I, i - 1, j + 1, k + 2);
	else if (internal_cell(I, i + 1, j - 1, k + 2))
		return POR_IND(I, i + 1, j - 1, k + 2);
	else if (internal_cell(I, i - 2, j - 1, k + 1))
		return POR_IND(I, i - 2, j - 1, k + 1);
	else if (internal_cell(I, i - 2, j + 1, k - 1))
		return POR_IND(I, i - 2, j + 1, k - 1);
	else if (internal_cell(I, i - 1, j - 2, k + 1))
		return POR_IND(I, i - 1, j - 2, k + 1);
	else if (internal_cell(I, i + 1, j - 2, k - 1))
		return POR_IND(I, i + 1, j - 2, k - 1);
	else if (internal_cell(I, i - 1, j + 1, k - 2))
		return POR_IND(I, i - 1, j + 1, k - 2);
	else if (internal_cell(I, i + 1, j - 1, k - 2))
		return POR_IND(I, i + 1, j - 1, k - 2);
	else if (internal_cell(I, i + 2, j - 1, k - 1))
		return POR_IND(I, i + 2, j - 1, k - 1);
	else if (internal_cell(I, i - 1, j + 2, k - 1))
		return POR_IND(I, i - 1, j + 2, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k + 2))
		return POR_IND(I, i - 1, j - 1, k + 2);
	else if (internal_cell(I, i - 2, j - 1, k - 1))
		return POR_IND(I, i - 2, j - 1, k - 1);
	else if (internal_cell(I, i - 1, j - 2, k - 1))
		return POR_IND(I, i - 1, j - 2, k - 1);
	else if (internal_cell(I, i - 1, j - 1, k - 2))
		return POR_IND(I, i - 1, j - 1, k - 2);
	else if (internal_cell(I, i + 1, j + 2, k + 2))
		return POR_IND(I, i + 1, j + 2, k + 2);
	else if (internal_cell(I, i + 2, j + 1, k + 2))
		return POR_IND(I, i + 2, j + 1, k + 2);
	else if (internal_cell(I, i + 2, j + 2, k + 1))
		return POR_IND(I, i + 2, j + 2, k + 1);
	else if (internal_cell(I, i - 1, j + 2, k + 2))
		return POR_IND(I, i - 1, j + 2, k + 2);
	else if (internal_cell(I, i + 2, j - 1, k + 2))
		return POR_IND(I, i + 2, j - 1, k + 2);
	else if (internal_cell(I, i + 2, j + 2, k - 1))
		return POR_IND(I, i + 2, j + 2, k - 1);
	else if (internal_cell(I, i + 1, j - 2, k + 2))
		return POR_IND(I, i + 1, j - 2, k + 2);
	else if (internal_cell(I, i + 1, j + 2, k - 2))
		return POR_IND(I, i + 1, j + 2, k - 2);
	else if (internal_cell(I, i - 2, j + 1, k + 2))
		return POR_IND(I, i - 2, j + 1, k + 2);
	else if (internal_cell(I, i + 2, j + 1, k - 2))
		return POR_IND(I, i + 2, j + 1, k - 2);
	else if (internal_cell(I, i - 2, j + 2, k + 1))
		return POR_IND(I, i - 2, j + 2, k + 1);
	else if (internal_cell(I, i + 2, j - 2, k + 1))
		return POR_IND(I, i + 2, j - 2, k + 1);
	else if (internal_cell(I, i - 1, j - 2, k + 2))
		return POR_IND(I, i - 1, j - 2, k + 2);
	else if (internal_cell(I, i - 1, j + 2, k - 2))
		return POR_IND(I, i - 1, j + 2, k - 2);
	else if (internal_cell(I, i - 2, j - 1, k + 2))
		return POR_IND(I, i - 2, j - 1, k + 2);
	else if (internal_cell(I, i + 2, j - 1, k - 2))
		return POR_IND(I, i + 2, j - 1, k - 2);
	else if (internal_cell(I, i - 2, j + 2, k - 1))
		return POR_IND(I, i - 2, j + 2, k - 1);
	else if (internal_cell(I, i + 2, j - 2, k - 1))
		return POR_IND(I, i + 2, j - 2, k - 1);
	else if (internal_cell(I, i + 1, j - 2, k - 2))
		return POR_IND(I, i + 1, j - 2, k - 2);
	else if (internal_cell(I, i - 2, j + 1, k - 2))
		return POR_IND(I, i - 2, j + 1, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k + 1))
		return POR_IND(I, i - 2, j - 2, k + 1);
	else if (internal_cell(I, i - 1, j - 2, k - 2))
		return POR_IND(I, i - 1, j - 2, k - 2);
	else if (internal_cell(I, i - 2, j - 1, k - 2))
		return POR_IND(I, i - 2, j - 1, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k - 1))
		return POR_IND(I, i - 2, j - 2, k - 1);
	else if (internal_cell(I, i + 2, j + 2, k + 2))
		return POR_IND(I, i + 2, j + 2, k + 2);
	else if (internal_cell(I, i - 2, j + 2, k + 2))
		return POR_IND(I, i - 2, j + 2, k + 2);
	else if (internal_cell(I, i + 2, j - 2, k + 2))
		return POR_IND(I, i + 2, j - 2, k + 2);
	else if (internal_cell(I, i + 2, j + 2, k - 2))
		return POR_IND(I, i + 2, j + 2, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k + 2))
		return POR_IND(I, i - 2, j - 2, k + 2);
	else if (internal_cell(I, i - 2, j + 2, k - 2))
		return POR_IND(I, i - 2, j + 2, k - 2);
	else if (internal_cell(I, i + 2, j - 2, k - 2))
		return POR_IND(I, i + 2, j - 2, k - 2);
	else if (internal_cell(I, i - 2, j - 2, k - 2))
		return POR_IND(I, i - 2, j - 2, k - 2);
}

int AREA_IND(in *I, int i, int j, int k, int s)
{
	return (6 * I->ind_cell_multipl[i * I->ny + j] + s);
}

int NORMAL_IND(in *I, int p, int i, int j, int k, int s)
{
	return (18 * I->ind_cell_multipl[i * I->ny + j] + s * 3 + p);
}

int VOLUME_IND(in *I, int i, int j, int k)
{
	return (I->ind_cell_multipl[i * I->ny + j]);
}

int check_for_corrupt_cell(in *I, int i, int j, int k)
{
	if ((i >= I->nx) || (i < 0) ||
		(j >= I->ny) || (j < 0) ||
		(k >= I->nz) || (k < 0) ||
		(I->ind_cell_multipl[i * I->ny + j] == -1)) {
			printf("!!!!!!Operating with corrupt cell!!!!!!\n");
			return 1;
	}
	return 0;
}

int well(in *I, int i, int j, int k)
{
	if ((injection_well(I, i, j, k)) || (production_well(I, i, j, k)))
		return 1;
	else
		return 0;
}

int production_well(in *I, int i, int j, int k)
{
	if ((i + I->ind_start_region_proc[0] == 0) && (j + I->ind_start_region_proc[1] == 0))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == 0) && (j + I->ind_start_region_proc[1] == I->gl_ny - 1))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx - 1) && (j + I->ind_start_region_proc[1] == 0))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx - 1) && (j + I->ind_start_region_proc[1] == I->gl_ny - 1))
		return 1;
	else
		return 0;
}

int injection_well(in *I, int i, int j, int k)
{
	if ((i + I->ind_start_region_proc[0] == I->gl_nx / 2) && (j + I->ind_start_region_proc[1] == I->gl_ny / 2))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx / 2 - 1) && (j + I->ind_start_region_proc[1] == I->gl_ny / 2))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx / 2) && (j + I->ind_start_region_proc[1] == I->gl_ny / 2 - 1))
		return 1;
	else if ((i + I->ind_start_region_proc[0] == I->gl_nx / 2 - 1) && (j + I->ind_start_region_proc[1] == I->gl_ny / 2 - 1))
		return 1;
	else
		return 0;
}

int write_to_A_csr(in *I, int p_eq, int i_eq, int j_eq, int k_eq, int p, int i, int j, int k, int s, double value)
{
	if (value != value)
		printf("Writing %lf to matrix A in (%d, %d, %d, %d, %d, %d, %d, %d)\n", value, p_eq, i_eq, j_eq, k_eq, p, i, j, k);
	if ((i >= 0) && (i < I->nx) && (j >= 0) && (j < I->ny) && (k >= 0) && (k < I->nz) && (I->ind_cell_multipl[i * I->ny + j] != -1) && (A_IND_S_SWITCH(I, i, j, k, s))) {
		int l, m, column, row;
		if (!((i_eq >= 0) && (i_eq < I->nx) && (j_eq >= 0) && (j_eq < I->ny) && (k_eq >= 0) && (k_eq < I->nz) && (I->ind_cell_multipl[i_eq * I->ny + j_eq] != -1)))
			printf("something very strange in ./src/lib/utils.c:87\n");
		if (s == -1) {
			column = A_IND(I, p, i, j, k);
		} else {
			column = A_IND_S(I, p, i, j, k, s);
		}
		row = A_IND(I, p_eq, i_eq, j_eq, k_eq);
		if (I->flag_first_time_step) {
			for (l = I->Aiptr_csr[row]; l <= I->A_ind_current; l++) {
				if (I->Ajptr_csr[l] == -1) {
					I->Aelem_csr[l] = value;
					I->Ajptr_csr[l] = column;
					break;
				} else if (I->Ajptr_csr[l] == column) {
					I->Aelem_csr[l] += value;
					break;
				} else if (I->Ajptr_csr[l] > column) {
					I->A_ind_current++;
					if (I->A_ind_current == I->non_zero_elem) {
						I->non_zero_elem++;
						if ((I->Aelem_csr = (double *) realloc(I->Aelem_csr, I->non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((I->Ajptr_csr = (int *) realloc(I->Ajptr_csr, I->non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					for (m = I->A_ind_current; m > l; m--) {
						I->Aelem_csr[m] = I->Aelem_csr[m - 1];
						I->Ajptr_csr[m] = I->Ajptr_csr[m - 1];
					}
					I->Aelem_csr[l] = value;
					I->Ajptr_csr[l] = column;
					break;
				} else if (l == I->A_ind_current) {
					I->A_ind_current++;
					if (I->A_ind_current == I->non_zero_elem) {
						I->non_zero_elem++;
						if ((I->Aelem_csr = (double *) realloc(I->Aelem_csr, I->non_zero_elem * sizeof(double))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
						if ((I->Ajptr_csr = (int *) realloc(I->Ajptr_csr, I->non_zero_elem * sizeof(int))) == NULL) {
							printf("Memory error\n");
							return 1;
						}
					}
					I->Aelem_csr[I->A_ind_current] = value;
					I->Ajptr_csr[I->A_ind_current] = column;
					break;
				}
			}
		} else {
			for (l = I->Aiptr_csr[row]; l < I->Aiptr_csr[row + 1]; l++) {
				if (I->Ajptr_csr[l] == column) {
					I->Aelem_csr[l] += value;
					break;
				}
			}
		}
	} else {
		I->B[A_IND(I, p_eq, i_eq, j_eq, k_eq)] -= value * I->B_prev[B_IND(I, p, i, j, k)];
	}
	return 0;
}

int A_IND_S_SWITCH(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			if ((i + 1 >= I->nx) || (I->ind_cell_multipl[(i + 1) * I->ny + j] == -1))
				return 0;
			else
				return 1;
		case 1:
			if ((i - 1 < 0) || (I->ind_cell_multipl[(i - 1) * I->ny + j] == -1))
				return 0;
			else
				return 1;
		case 2:
			if ((j + 1 >= I->ny) || (I->ind_cell_multipl[i * I->ny + j + 1] == -1))
				return 0;
			else
				return 1;
		case 3:
			if ((j - 1 < 0) || (I->ind_cell_multipl[i * I->ny + j - 1] == -1))
				return 0;
			else
				return 1;
		case 4:
			if (k + 1 >= I->nz)
				return 0;
			else
				return 1;
		case 5:
			if (k - 1 < 0)
				return 0;
			else
				return 1;
	}
	return 1;
}

int A_IND_S(in *I, int p, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return A_IND(I, p, i + 1, j, k);
		case 1:
			return A_IND(I, p, i - 1, j, k);
		case 2:
			return A_IND(I, p, i, j + 1, k);
		case 3:
			return A_IND(I, p, i, j - 1, k);
		case 4:
			return A_IND(I, p, i, j, k + 1);
		case 5:
			return A_IND(I, p, i, j, k - 1);
	}
}

int cell_of_computation_domain(in *I, int i, int j, int k)
{
	if ((i >= - I->stencil_size) && (i < I->nx + I->stencil_size) &&
		(j >= - I->stencil_size) && (j < I->ny + I->stencil_size) &&
		(k >= - I->stencil_size) && (k < I->nz + I->stencil_size))
			return 1;
	else
			return 0;
}

int internal_cell(in *I, int i, int j, int k)
{
	if ((i >= 0) && (i < I->nx) && (j >= 0) && (j < I->ny) && (I->ind_cell_multipl[i * I->ny + j] != -1) && (k >= 0) && (k < I->nz))
		return 1;
	else
		return 0;

}

int boundary_cell(in *I, int i, int j, int k)
{
	if (cell_of_computation_domain(I, i, j, k) &&
		(I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
		(!(internal_cell(I, i, j, k))))
			return 1;
	else
			return 0;

}

int count_neighbor_internal_cells(in *I, int i, int j, int k)
{
	int x = 0;
	if (boundary_cell(I, i - 1, j, k)) x++;
	if (boundary_cell(I, i + 1, j, k)) x++;
	if (boundary_cell(I, i, j - 1, k)) x++;
	if (boundary_cell(I, i, j + 1, k)) x++;
	return x;
}

int count_second_order_neighbor_internal_cells(in *I, int i, int j, int k)
{
	int x = 0;
	if (boundary_cell(I, i - 2, j, k)) x++;
	if (boundary_cell(I, i + 2, j, k)) x++;
	if (boundary_cell(I, i, j - 2, k)) x++;
	if (boundary_cell(I, i, j + 2, k)) x++;
	return x;
}

int count_other_corner_neighbor_internal_cells(in *I, int i, int j, int k)
{
	int x = 0;
	if (boundary_cell(I, i - 1, j - 1, k)) x++;
	if (boundary_cell(I, i - 1, j + 1, k)) x++;
	if (boundary_cell(I, i + 1, j - 1, k)) x++;
	if (boundary_cell(I, i + 1, j + 1, k)) x++;
	if (boundary_cell(I, i - 1, j - 2, k)) x++;
	if (boundary_cell(I, i - 1, j + 2, k)) x++;
	if (boundary_cell(I, i + 1, j - 2, k)) x++;
	if (boundary_cell(I, i + 1, j + 2, k)) x++;
	if (boundary_cell(I, i - 2, j - 1, k)) x++;
	if (boundary_cell(I, i - 2, j + 1, k)) x++;
	if (boundary_cell(I, i + 2, j - 1, k)) x++;
	if (boundary_cell(I, i + 2, j + 1, k)) x++;
	if (boundary_cell(I, i - 2, j - 2, k)) x++;
	if (boundary_cell(I, i - 2, j + 2, k)) x++;
	if (boundary_cell(I, i + 2, j - 2, k)) x++;
	if (boundary_cell(I, i + 2, j + 2, k)) x++;
	return x;
}

double max2(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
}

double min2(double a, double b)
{
	if (a < b)
		return a;
	else
		return b;
}

double max3(double a, double b, double c)
{
	if ((a > b) && (a > c))
		return a;
	else if (b > c)
		return b;
	else
		return c;
}

double min3(double a, double b, double c)
{
	if ((a < b) && (a < c))
		return a;
	else if (b < c)
		return b;
	else
		return c;
}

int nan_check(in *I, double value, int p, int pr, int i, int j, int k, const char func_name[100])
{
	if (isnan(value) || fabs(isinf(value))) {
		I->nan_flag = 1;
		printf("%s (", func_name);
		if (p != -1)
			printf("p = %d, ", p);
		if (pr != -1)
			printf("pr = %d, ", pr);
		printf("i = %d, j = %d, k = %d) = %lf\n", i, j, k, value);
		return 1;
	} else
		return 0;
}

int negative_or_zero_num_check(in *I, double value, int p, int pr, int i, int j, int k, const char func_name[100])
{
	if (value <= 0) {
		I->negative_num_flag = 1;
		printf("%s (", func_name);
		if (p != -1)
			printf("p = %d, ", p);
		if (pr != -1)
			printf("pr = %d, ", pr);
		printf("i = %d, j = %d, k = %d) = %lf\n", i, j, k, value);
		return 1;
	} else
		return 0;
}

#if AVALANCHE
double density(in *I, int i, int j, int k)
{
	//return phase_fraction(I, i, j, k) * I->density_snow + (1 - phase_fraction(I, i, j, k)) * I->density_air;
	return phase_fraction(I, i, j, k);
}

double velocity(in *I, int p, int i, int j, int k)
{
	return I->B_prev[B_IND(I, p, i, j, k)];
}

double phase_fraction(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 3, i, j, k)];
}

double pressure(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 4, i, j, k)];
}

double velocity_on_face(in *I, int p, int i, int j, int k, int s)
{
	int a, b;
	switch (s) {
		case 0:
			return (velocity(I, p, i, j, k) + velocity(I, p, i + 1, j, k)) / 2;
		case 1:
			return (velocity(I, p, i, j, k) + velocity(I, p, i - 1, j, k)) / 2;
		case 2:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j + 1, k)) / 2;
		case 3:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j - 1, k)) / 2;
		case 4:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j, k + 1)) / 2;
		case 5:
			return (velocity(I, p, i, j, k) + velocity(I, p, i, j, k - 1)) / 2;
	}
}

double density_on_face(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (density(I, i, j, k) + density(I, i + 1, j, k)) / 2;
		case 1:
			return (density(I, i, j, k) + density(I, i - 1, j, k)) / 2;
		case 2:
			return (density(I, i, j, k) + density(I, i, j + 1, k)) / 2;
		case 3:
			return (density(I, i, j, k) + density(I, i, j - 1, k)) / 2;
		case 4:
			return (density(I, i, j, k) + density(I, i, j, k + 1)) / 2;
		case 5:
			return (density(I, i, j, k) + density(I, i, j, k - 1)) / 2;
	}
}

double strain_rate_on_face(in *I, int i, int j, int k, int s, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity_on_face(I, 0, i + 1, j, k, s) - velocity_on_face(I, 0, i - 1, j, k, s)) / (2 * I->dx[0]);
		case 1:
			return 0.5 * ((velocity_on_face(I, 0, i, j + 1, k, s) - velocity_on_face(I, 0, i, j - 1, k, s)) / (2 * I->dx[1]) +
					(velocity_on_face(I, 1, i + 1, j, k, s) - velocity_on_face(I, 1, i - 1, j, k, s)) / (2 * I->dx[0]));
		case 2:
			return 0.5 * ((velocity_on_face(I, 0, i, j, k + 1, s) - velocity_on_face(I, 0, i, j, k - 1, s)) / (2 * I->dx[2]) +
					(velocity_on_face(I, 2, i + 1, j, k, s) - velocity_on_face(I, 2, i - 1, j, k, s)) / (2 * I->dx[0]));
		case 3:
			return (velocity_on_face(I, 1, i, j + 1, k, s) - velocity_on_face(I, 1, i, j - 1, k, s)) / (2 * I->dx[1]);
		case 5:
			return 0.5 * ((velocity_on_face(I, 1, i, j, k + 1, s) - velocity_on_face(I, 1, i, j, k - 1, s)) / (2 * I->dx[2]) +
					(velocity_on_face(I, 2, i, j + 1, k, s) - velocity_on_face(I, 2, i, j - 1, k, s)) / (2 * I->dx[1]));
		case 8:
			return (velocity_on_face(I, 2, i, j, k + 1, s) - velocity_on_face(I, 2, i, j, k - 1, s)) / (2 * I->dx[2]);
	}
}

double shear_rate_on_face(in *I, int i, int j, int k, int s)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate_on_face(I, i, j, k, s, m, n) * strain_rate_on_face(I, i, j, k, s, m, n);
		}
	}
	return x;
}

double strain_rate(in *I, int i, int j, int k, int m, int n)
{
	switch (m + n + m * n) {
		case 0:
			return (velocity(I, 0, i + 1, j, k) - velocity(I, 0, i - 1, j, k)) / (2 * I->dx[0]);
		case 1:
			return 0.5 * ((velocity(I, 0, i, j + 1, k) - velocity(I, 0, i, j - 1, k)) / (2 * I->dx[1]) +
				(velocity(I, 1, i + 1, j, k) - velocity(I, 1, i - 1, j, k)) / (2 * I->dx[0]));
		case 2:
			return 0.5 * ((velocity(I, 0, i, j, k + 1) - velocity(I, 0, i, j, k - 1)) / (2 * I->dx[2]) +
				(velocity(I, 2, i + 1, j, k) - velocity(I, 2, i - 1, j, k)) / (2 * I->dx[0]));
		case 3:
			return (velocity(I, 1, i, j + 1, k) - velocity(I, 1, i, j - 1, k)) / (2 * I->dx[1]);
		case 5:
			return 0.5 * ((velocity(I, 1, i, j, k + 1) - velocity(I, 1, i, j, k - 1)) / (2 * I->dx[2]) +
				(velocity(I, 2, i, j + 1, k) - velocity(I, 2, i, j - 1, k)) / (2 * I->dx[1]));
		case 8:
			return (velocity(I, 2, i, j, k + 1) - velocity(I, 2, i, j, k - 1)) / (2 * I->dx[2]);
	}
}

double shear_rate(in *I, int i, int j, int k)
{
	int m, n;
	double x = 0;
	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n ++) {
			x += strain_rate(I, i, j, k, m, n) * strain_rate(I, i, j, k, m, n);
		}
	}
	return x;
}

double phase_fraction_on_face(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i + 1, j, k)) / 2;
		case 1:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i - 1, j, k)) / 2;
		case 2:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j + 1, k)) / 2;
		case 3:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j - 1, k)) / 2;
		case 4:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j, k + 1)) / 2;
		case 5:
			return (phase_fraction(I, i, j, k) + phase_fraction(I, i, j, k - 1)) / 2;
	}
}

double effective_viscosity_on_face(in *I, int i, int j, int k, int s)
{
	double phase_fraction_on_face_tmp = phase_fraction_on_face(I, i, j, k, s);
	double shear_rate_on_face_tmp = shear_rate_on_face(I, i, j, k, s);
	double x = I->k_viscosity_air * (1 - phase_fraction_on_face_tmp);
	if (shear_rate_on_face_tmp <= I->shear_rate_0) {
		x += I->limiting_viscosity_snow * phase_fraction_on_face_tmp;
	} else {
		x += (I->k_viscosity_snow * pow(shear_rate_on_face_tmp, I->flow_index - 1) + I->yield_stress / shear_rate_on_face_tmp) * phase_fraction_on_face_tmp;
	}
	return x;
}

double pressure_on_face(in *I, int i, int j, int k, int s)
{
	switch (s) {
		case 0:
			return (pressure(I, i, j, k) + pressure(I, i + 1, j, k)) / 2;
		case 1:
			return (pressure(I, i, j, k) + pressure(I, i - 1, j, k)) / 2;
		case 2:
			return (pressure(I, i, j, k) + pressure(I, i, j + 1, k)) / 2;
		case 3:
			return (pressure(I, i, j, k) + pressure(I, i, j - 1, k)) / 2;
		case 4:
			return (pressure(I, i, j, k) + pressure(I, i, j, k + 1)) / 2;
		case 5:
			return (pressure(I, i, j, k) + pressure(I, i, j, k - 1)) / 2;
	}
}

int check_conservation_of_mass(in *I)
{
	int i, j, k;
	double mass = 0;
	for (k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if ((I->ind_boundary_cells[(i + I->stencil_size) * (I->ny + 2 * I->stencil_size) + j + I->stencil_size] != -1) &&
					((!((i >= 0) && (i < I->nx) &&
					   (j >= 0) && (j < I->ny) && 
					   (I->ind_cell_multipl[i * I->ny + j] != -1))) ||
					 (k < 0) || (k >= I->nz))) {
							I->mass_quantity -= phase_fraction(I, i, j, k);
				}
			}
		}
	}
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					mass += phase_fraction(I, i, j, k);
				}
			}
		}
	}
	if (mass == I->mass_quantity) {
		return 0;
	} else {
		return 1;
	}
}

int barotropy_pressure(in *I)
{
	int i, j, k;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					if (check_for_corrupt_cell(I, i, j, k)) return 1;
					I->B[A_IND(I, 4, i, j, k)] = I->pressure_atmosphere * (
						1 +
						0.001 * (I->B[A_IND(I, 3, i, j, k)] - I->density_snow) / I->density_snow +
						0 * (I->B[A_IND(I, 3, i, j, k)] - I->density_snow) * (I->B[A_IND(I, 3, i, j, k)] - I->density_snow) / (I->density_snow * I->density_snow));
				}
			}
		}
	}
	return 0;
}

int barotropy_density(in *I)
{
	int i, j, k;
	double c1 = 1, c2 = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					if (check_for_corrupt_cell(I, i, j, k)) return 1;
					I->B[A_IND(I, 3, i, j, k)] = (I->density_snow / I->pressure_atmosphere) * (
						I->pressure_atmosphere +
						c1 * (I->B[A_IND(I, 4, i, j, k)] - I->pressure_atmosphere) +
						c2 * 0.5 * (I->B[A_IND(I, 4, i, j, k)] - I->pressure_atmosphere) *
						(I->B[A_IND(I, 4, i, j, k)] - I->pressure_atmosphere));
				}
			}
		}
	}
	return 0;
}
#endif

#if TERMOGAS
double saturation(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2))) {
		printf("Error saturation index\n");
		return 0;
	}
	return I->B_prev[B_IND(I, p + 5, i, j, k)] + I->residual_saturation[p];
}

double concentration(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2) || (p == 3))) {
		printf("Error concentration index\n");
		return 0;
	}
	return I->B_prev[B_IND(I, p, i, j, k)];
}

double pressure(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 4, i, j, k)];
}

double temperature_flow(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 8, i, j, k)] + I->initial_temperature;
}

double temperature_environment(in *I, int i, int j, int k)
{
	return I->B_prev[B_IND(I, 9, i, j, k)] + I->initial_temperature;
}

double density_t_func(in *I, int p, int i, int j, int k)
{
	double tmp;
	if ((p == 0) || (p == 1))
		tmp = ((I->density_0[p] + pow(I->density_coef_a[p], -2) * (pressure(I, i, j, k) - I->pressure_0)) /
			(1 + I->density_coef_beta[p] * (temperature_flow(I, i, j, k) - I->temperature_0)));
	else if (p == 2)
		tmp = (pressure(I, i, j, k) * (concentration(I, 0, i, j, k) * I->molar_weight[0] +
										concentration(I, 1, i, j, k) * I->molar_weight[1] +
										concentration(I, 2, i, j, k) * I->molar_weight[2] +
										concentration(I, 3, i, j, k) * I->molar_weight[3]) / (I->R * temperature_flow(I, i, j, k)));
	else {
		printf("Error density index\n");
		return 0;
	}
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	negative_or_zero_num_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double density_t(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2))
		return I->array_of_parameters[PARAM_IND(I, p, i, j, k)];
	else {
		printf("Error density index\n");
		return 0;
	}
}

double two_phase_relative_permeability_func(in *I, int p, int pr, int i, int j, int k)
{
	double S, tmp;
	if ((p == 0) && (pr == 1)) {
		S = (saturation(I, 0, i, j, k) - I->residual_saturation[0]);// / (1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2]);
		tmp = (pow(S, 4));
	}
	else if ((p == 1) && (pr == 0)) {
		S = (saturation(I, 0, i, j, k) - I->residual_saturation[0]);// / (1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2]);
		tmp = ((1 - S) * (1 - S) * (1 - S * S));
	}
	else if ((p == 0) && (pr == 2)) {
		S = (saturation(I, 2, i, j, k) - I->residual_saturation[2]);// / (1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2]);
		tmp = ((1 - S) * (1 - S) * (1 - S * S));
	}
	else if ((p == 2) && (pr == 0)) {
		S = (saturation(I, 2, i, j, k) - I->residual_saturation[2]);// / (1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2]);
		tmp = (pow(S, 4));
	}
	else if ((p == 1) && (pr == 2)) {
		S = (saturation(I, 1, i, j, k) - I->residual_saturation[1]);// / (1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2]);
		tmp = (pow(S, 4));
	}
	else if ((p == 2) && (pr == 1)) {
		S = (saturation(I, 1, i, j, k) - I->residual_saturation[1]);// / (1 - I->residual_saturation[0] - I->residual_saturation[1] - I->residual_saturation[2]);
		tmp = ((1 - S) * (1 - S) * (1 - S * S));
	} else {
		printf("Error two phase relative permeability index\n");
		return 0;
	}
	nan_check(I, tmp, p, pr, i, j, k, __func__);
	return tmp;
}

double two_phase_relative_permeability(in *I, int p, int pr, int i, int j, int k)
{
	double S;
	if ((p == 0) && (pr == 1))
		return I->array_of_parameters[PARAM_IND(I, 3, i, j, k)];
	else if ((p == 1) && (pr == 0))
		return I->array_of_parameters[PARAM_IND(I, 4, i, j, k)];
	else if ((p == 0) && (pr == 2))
		return I->array_of_parameters[PARAM_IND(I, 5, i, j, k)];
	else if ((p == 2) && (pr == 0))
		return I->array_of_parameters[PARAM_IND(I, 6, i, j, k)];
	else if ((p == 1) && (pr == 2))
		return I->array_of_parameters[PARAM_IND(I, 7, i, j, k)];
	else if ((p == 2) && (pr == 1))
		return I->array_of_parameters[PARAM_IND(I, 8, i, j, k)];
	else {
		printf("Error two phase relative permeability index\n");
		return 0;
	}
}

double relative_permeability_func(in *I, int p, int i, int j, int k)
{
	double tmp;
	if (p == 0) {
		if (saturation(I, 1, i, j, k) - I->residual_saturation[1] + saturation(I, 2, i, j, k) - I->residual_saturation[2] == 0)
			tmp = 0;
		else
			tmp =  ((saturation(I, 1, i, j, k) - I->residual_saturation[1]) * two_phase_relative_permeability(I, 0, 1, i, j, k) +
					(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * two_phase_relative_permeability(I, 0, 2, i, j, k)) /
					(saturation(I, 1, i, j, k) - I->residual_saturation[1] + saturation(I, 2, i, j, k) - I->residual_saturation[2]);
		//tmp =  ((saturation(I, 1, i, j, k) - I->residual_saturation[1] + I->epsilon) * pow(two_phase_relative_permeability(I, 0, 1, i, j, k), 0.25) +
		//		(saturation(I, 2, i, j, k) - I->residual_saturation[2] + I->epsilon) * pow(two_phase_relative_permeability(I, 0, 2, i, j, k), 0.25)) /
		//		(saturation(I, 1, i, j, k) - I->residual_saturation[1] + saturation(I, 2, i, j, k) - I->residual_saturation[2] + I->epsilon);
	} else if (p == 1) {
		if (saturation(I, 0, i, j, k) - I->residual_saturation[0] + saturation(I, 2, i, j, k) - I->residual_saturation[2] == 0)
			tmp = 0;
		else
			tmp =  ((saturation(I, 0, i, j, k) - I->residual_saturation[0]) * two_phase_relative_permeability(I, 1, 0, i, j, k) +
					(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * two_phase_relative_permeability(I, 1, 2, i, j, k)) /
					(saturation(I, 0, i, j, k) - I->residual_saturation[0] + saturation(I, 2, i, j, k) - I->residual_saturation[2]);
	} else if (p == 2) {
		if (saturation(I, 0, i, j, k) - I->residual_saturation[0] + saturation(I, 1, i, j, k) - I->residual_saturation[1] == 0)
			tmp = 0;
		else
			tmp =  ((saturation(I, 0, i, j, k) - I->residual_saturation[0]) * two_phase_relative_permeability(I, 2, 0, i, j, k) +
					(saturation(I, 1, i, j, k) - I->residual_saturation[1]) * two_phase_relative_permeability(I, 2, 1, i, j, k)) /
					(saturation(I, 0, i, j, k) - I->residual_saturation[0] + saturation(I, 1, i, j, k) - I->residual_saturation[1]);
	} else {
		printf("Error relative permeability index\n");
		return 0;
	}
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double relative_permeability(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2))
		return I->array_of_parameters[PARAM_IND(I, 9 + p, i, j, k)];
	else {
		printf("Error relative permeability index\n");
		return 0;
	}
}

double viscosity_gas_component_func(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2) || (p == 3))) {
		printf("Error viscosity gas index\n");
		return 0;
	}
	double tmp =  I->viscosity_coef_A_gas[p] *
		(I->temperature_0_gas[p] + I->viscosity_coef_C_gas[p]) *
		pow(temperature_flow(I, i, j, k) / I->temperature_0_gas[p], 3 / 2) /
		(temperature_flow(I, i, j, k) + I->viscosity_coef_C_gas[p]);
	negative_or_zero_num_check(I, tmp, p, -1, i, j, k, __func__);
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double viscosity_gas_component(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2) || (p == 3))
		return I->array_of_parameters[PARAM_IND(I, 12 + p, i, j, k)];
	else {
		printf("Error viscosity gas index\n");
		return 0;
	}
}

double molar_fraction_func(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 1) || (p == 2) || (p == 3))) {
		printf("Error molar fraction index\n");
		return 0;
	}
	double x = 0;
	int l;
	for (l = 0; l < 4; l++)
		x += concentration(I, l, i, j, k) / I->molar_weight[l];
	double tmp = concentration(I, p, i, j, k) / (I->molar_weight[p] * x);
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double molar_fraction(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2) || (p == 3))
		return I->array_of_parameters[PARAM_IND(I, 16 + p, i, j, k)];
	else {
		printf("Error molar fraction index\n");
		return 0;
	}
}

double viscosity_water(in *I, int i, int j, int k)
{
	double tmp;
	tmp = 0.001 / (0.14 + (temperature_flow(I, i, j, k) - 273.15) / 30.0 + 0.000009 * (temperature_flow(I, i, j, k) - 273.15) * (temperature_flow(I, i, j, k) - 273.15));
	negative_or_zero_num_check(I, tmp, -1, -1, i, j, k, __func__);
	if (temperature_flow(I, i, j, k) - 273.15 < 0)
		printf("WTF?! Temperature flow is less then 273.15 K!!!\n");
	return tmp;
}

double viscosity_oil(in *I, int i, int j, int k)
{
	double tmp;
	tmp = density_t(I, 1, i, j, k) * 0.000001 * (pow((1000.0 + 0.6), pow((30.0 / (temperature_flow(I, i, j, k) - 273.15)), 4.0)) - 0.6);
	negative_or_zero_num_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double viscosity_gas(in *I, int i, int j, int k)
{
	double tmp = 0;
	for (int p = 0; p < 4; p++)
		tmp += viscosity_gas_component(I, p, i, j, k) * molar_fraction(I, p, i, j, k);
	negative_or_zero_num_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double viscosity_func(in *I, int p, int i, int j, int k)
{
	double tmp = 0;
	if (p == 0)
		tmp = viscosity_water(I, i, j, k);
	else if (p == 1)
		tmp = viscosity_oil(I, i, j, k);
	else if (p == 2)
		tmp = viscosity_gas(I, i, j, k);
	else
		printf("Error viscosity index\n");
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double viscosity(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2))
		return I->array_of_parameters[PARAM_IND(I, 20 + p, i, j, k)];
	else {
		printf("Error viscosity index\n");
		return 0;
	}
}

double Darsi_M_coef_phases_func(in *I, int p, int i, int j, int k)
{
	double tmp = 0;
	if ((p == 0) || (p == 1) || (p == 2))
		tmp = density_t(I, p, i, j, k) * I->permeability * relative_permeability(I, p, i, j, k) / viscosity(I, p, i, j, k);
	else
		printf("Error Darsi M coefficient phases index\n");
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double Darsi_M_coef_phases(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2))
		return I->array_of_parameters[PARAM_IND(I, 23 + p, i, j, k)];
	else {
		printf("Error Darsi M coefficient phases index\n");
		return 0;
	}
}

double Darsi_M_coef_func(in *I, int i, int j, int k)
{
	double tmp = 0;
	for (int p = 0; p < 3; p++)
		tmp += Darsi_M_coef_phases(I, p, i, j, k);
	if (relative_permeability(I, 0, i, j, k) + relative_permeability(I, 1, i, j, k) + relative_permeability(I, 2, i, j, k) == 0) {
		printf("Relative permeability sum is equal to zero!!!\n");
		printf("saturation(I, 0, %d, %d, %d) = %lf, saturation(I, 1, %d, %d, %d) = %lf, saturation(I, 2, %d, %d, %d) = %lf.\n",\
			i, j, k, saturation(I, 0, i, j, k), i, j, k, saturation(I, 1, i, j, k), i, j, k, saturation(I, 2, i, j, k));
	}
	nan_check(I, tmp, -1, -1, i, j, k, __func__);
	negative_or_zero_num_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double Darsi_M_coef(in *I, int i, int j, int k)
{
	return I->array_of_parameters[PARAM_IND(I, 26, i, j, k)];
}

double capillary_pressure_derivative_by_saturation_func(in *I, int p, int i, int j, int k)
{
	if (!((p == 0) || (p == 2))) {
		printf("Error capillary pressure index\n");
		return 0;
	}
	double tmp = - I->capillary_pressure_at_maximum_saturation[p] *
		pow(1 - I->residual_saturation_two_phase[p], I->capillary_pressure_coef) *
		I->capillary_pressure_coef *
		pow(saturation(I, p, i, j, k), - I->capillary_pressure_coef - 1);
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double capillary_pressure_derivative_by_saturation(in *I, int p, int i, int j, int k)
{
	if (p == 0)
		return I->array_of_parameters[PARAM_IND(I, 27, i, j, k)];
	else if (p == 2)
		return I->array_of_parameters[PARAM_IND(I, 28, i, j, k)];
	else {
		printf("Error capillary pressure index\n");
		return 0;
	}
}

double grad_pressure_func(in *I, int pr, int i, int j, int k)
{
	if (!cell_of_computation_domain(I, i, j, k)) {
		printf("%s: Operating with the non-existent cell (%d, %d, %d).\n", __func__, i, j, k);
		return 0;
	}
	int ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	double tmp = 0;
	if (!cell_of_computation_domain(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))
		tmp = (pressure(I, i, j, k) - pressure(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / I->dx[pr];
	else if (!cell_of_computation_domain(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))
		tmp = (pressure(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - pressure(I, i, j, k)) / I->dx[pr];
	else
		tmp = (pressure(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - pressure(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
	nan_check(I, tmp, -1, pr, i, j, k, __func__);
	return tmp;
}

double grad_pressure(in *I, int pr, int i, int j, int k)
{
	if (!cell_of_computation_domain(I, i, j, k)) {
		printf("%s: Operating with the non-existent cell (%d, %d, %d).\n", __func__, i, j, k);
		return 0;
	}
	if ((pr == 0) || (pr == 1) || (pr == 2))
		return I->array_of_parameters[PARAM_IND(I, 29 + pr, i, j, k)];
	else {
		printf("Error gradient pressure index\n");
		return 0;
	}
}

double grad_saturation_func(in *I, int p, int pr, int i, int j, int k)
{
	if (!cell_of_computation_domain(I, i, j, k)) {
		printf("%s: Operating with the non-existent cell (%d, %d, %d).\n", __func__, i, j, k);
		return 0;
	}
	int ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	double tmp = 0;
	if (!cell_of_computation_domain(I, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]))
		tmp = (saturation(I, p, i, j, k) - saturation(I, p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / I->dx[pr];
	else if (!cell_of_computation_domain(I, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2]))
		tmp = (saturation(I, p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - saturation(I, p, i, j, k)) / I->dx[pr];
	else
		tmp = (saturation(I, p, i + ind_pr[0], j + ind_pr[1], k + ind_pr[2]) - saturation(I, p, i - ind_pr[0], j - ind_pr[1], k - ind_pr[2])) / (2 * I->dx[pr]);
	nan_check(I, tmp, p, pr, i, j, k, __func__);
	return tmp;
}

double grad_saturation(in *I, int p, int pr, int i, int j, int k)
{
	if (!cell_of_computation_domain(I, i, j, k)) {
		printf("%s: Operating with the non-existent cell (%d, %d, %d).\n", __func__, i, j, k);
		return 0;
	}
	if (((pr == 0) || (pr == 1) || (pr == 2)) && ((p == 0) || (p == 1) || (p == 2)))
		return I->array_of_parameters[PARAM_IND(I, 32 + p * 3 + pr, i, j, k)];
	else {
		printf("Error gradient pressure index\n");
		return 0;
	}
}

double coef_grad_saturation_func(in *I, int p, int pr, int i, int j, int k)
{
	double tmp = 0;
	if (p == 0) {
		if (pr == 0)
			tmp = Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k) - 1;
		else if (pr == 2)
			tmp = Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k);
		else
			printf("Error coef_grad_saturation index: pr = %d.\n", pr);
	} else if (p == 1) {
		if (pr == 0)
			tmp = Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k);
		else if (pr == 2)
			tmp = Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k);
		else
			printf("Error coef_grad_saturation index: pr = %d.\n", pr);
	} else if (p == 2) {
		if (pr == 0)
			tmp = Darsi_M_coef_phases(I, 0, i, j, k) / Darsi_M_coef(I, i, j, k);
		else if (pr == 2)
			tmp = Darsi_M_coef_phases(I, 2, i, j, k) / Darsi_M_coef(I, i, j, k) - 1;
		else
			printf("Error coef_grad_saturation index: pr = %d.\n", pr);
	} else
		printf("Error coef_grad_saturation index: p = %d.\n", p);
	nan_check(I, tmp, p, pr, i, j, k, __func__);
	return tmp;
}

double coef_grad_saturation(in *I, int p, int pr, int i, int j, int k)
{
	if (p == 0) {
		if (pr == 0)
			return I->array_of_parameters[PARAM_IND(I, 41, i, j, k)];
		else if (pr == 2)
			return I->array_of_parameters[PARAM_IND(I, 42, i, j, k)];
		else
			printf("Error coef_grad_saturation index.\n");
	} else if (p == 1) {
		if (pr == 0)
			return I->array_of_parameters[PARAM_IND(I, 43, i, j, k)];
		else if (pr == 2)
			return I->array_of_parameters[PARAM_IND(I, 44, i, j, k)];
		else
			printf("Error coef_grad_saturation index.\n");
	} else if (p == 2) {
		if (pr == 0)
			return I->array_of_parameters[PARAM_IND(I, 45, i, j, k)];
		else if (pr == 2)
			return I->array_of_parameters[PARAM_IND(I, 46, i, j, k)];
		else
			printf("Error coef_grad_saturation index.\n");
	} else
		printf("Error coef_grad_saturation index.\n");
}

double avarage_velocity_func(in *I, int p, int pr, int i, int j, int k) //p - oil, water, gas; pr - x1, x2, x3
{
	if (!((pr == 0) || (pr == 1) || (pr == 2))) {
		printf("Error avarage velocity index\n");
		return 0;
	}
	int ind_pr[3];
	ind_pr[0] = ind_pr[1] = ind_pr[2] = 0;
	ind_pr[pr] = 1;
	double tmp = - I->permeability * (relative_permeability(I, p, i, j, k) + I->epsilon) * grad_pressure(I, pr, i, j, k) / viscosity(I, p, i, j, k);
//	double tmp = - I->permeability * relative_permeability(I, p, i, j, k) *(
//		grad_pressure(I, pr, i, j, k) +
//		coef_grad_saturation(I, p, 0, i, j, k) * capillary_pressure_derivative_by_saturation(I, 0, i, j, k) * grad_saturation(I, 0, pr, i, j, k) +
//		coef_grad_saturation(I, p, 2, i, j, k) * capillary_pressure_derivative_by_saturation(I, 2, i, j, k) * grad_saturation(I, 2, pr, i, j, k)
//		) / viscosity(I, p, i, j, k);
	nan_check(I, tmp, p, pr, i, j, k, __func__);
	return tmp;
}

double avarage_velocity(in *I, int p, int pr, int i, int j, int k) //p - oil, water, gas; pr - x1, x2, x3
{
	if (((pr == 0) || (pr == 1) || (pr == 2)) && ((p == 0) || (p == 1) || (p == 2)))
		return I->array_of_parameters[PARAM_IND(I, 47 + p * 3 + pr, i, j, k)];
	else {
		printf("Error avarage velocity index\n");
		return 0;
	}
}

double rate_of_reaction_coef_func(in *I, int i, int j, int k)
{
	double tmp;
	if (temperature_flow(I, i, j, k) < I->threshold_temperature)
		tmp = 0;
	else
		tmp = I->frequency_factor * exp(- I->activation_temperature / temperature_flow(I, i, j, k));
	nan_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double rate_of_reaction_coef(in *I, int i, int j, int k)
{
	return I->array_of_parameters[PARAM_IND(I, 56, i, j, k)];
}

double rate_of_reaction_func(in *I, int i, int j, int k)
{
	double tmp = I->porousness[POR_IND(I, i, j, k)] * rate_of_reaction_coef(I, i, j, k) * (saturation(I, 1, i, j, k) - I->residual_saturation[1]) *
		(saturation(I, 2, i, j, k) - I->residual_saturation[2]) * concentration(I, 1, i, j, k) *
		pow(pressure(I, i, j, k) / I->pressure_activ, I->stoichiometric_coef_activ);
	nan_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double rate_of_reaction(in *I, int i, int j, int k)
{
	return I->array_of_parameters[PARAM_IND(I, 57, i, j, k)];
}

double mass_inflow_rate_func(in *I, int p, int i, int j, int k)
{
	double tmp;
	if (p == 0) // N2
		tmp = 0;
	else if (p == 1) // O2
		tmp = rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 2) // CO2
		tmp = rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 3) // H2O(g)
		tmp = rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p] - I->stoichiometric_coef_before[p]) * I->molar_weight[p];
	else if (p == 5) // water
		tmp = rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p - 1] - I->stoichiometric_coef_before[p - 1]) * I->molar_weight[3];
	else if (p == 6) // oil
		tmp = rate_of_reaction(I, i, j, k) * (I->stoichiometric_coef_after[p - 1] - I->stoichiometric_coef_before[p - 1]) * I->molar_weight[4];
	else if (p == 7) // gas
		tmp = mass_inflow_rate_func(I, 0, i, j, k) + mass_inflow_rate_func(I, 1, i, j, k) +
			mass_inflow_rate_func(I, 2, i, j, k) + mass_inflow_rate_func(I, 3, i, j, k);
	else {
		printf("Error mass inflow rate index\n");
		return 0;
	}
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double mass_inflow_rate(in *I, int p, int i, int j, int k)
{
	if ((p == 0) || (p == 1) || (p == 2) || (p == 3))
		return I->array_of_parameters[PARAM_IND(I, 58 + p, i, j, k)];
	else if ((p == 5) || (p == 6) || (p == 7))
		return I->array_of_parameters[PARAM_IND(I, 62 + p - 5, i, j, k)];
	else {
		printf("Error mass inflow rate index\n");
		return 0;
	}
}

double chemical_reaction_heat_flow_func(in *I, int i, int j, int k)
{
	int p;
	double tmp = 0;
	for (p = 0; p < 4; p++) {
		tmp -= mass_inflow_rate(I, p, i, j, k) * I->initial_enthalpy[p + 2];
	}
	for (p = 0; p < 2; p++) {
		tmp -= mass_inflow_rate(I, p + 5, i, j, k) * I->initial_enthalpy[p];
	}
	nan_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double chemical_reaction_heat_flow(in *I, int i, int j, int k)
{
	return I->array_of_parameters[PARAM_IND(I, 65, i, j, k)];
}

double density_derivative_by_pressure(in *I, int p, int i, int j, int k)
{
	double tmp;
	if ((p == 0) || (p == 1))
		tmp = pow(I->density_coef_a[p], -2) / (1 + I->density_coef_beta[p] * (temperature_flow(I, i, j, k) - I->temperature_0));
	else if (p == 2)
		tmp = (concentration(I, 0, i, j, k) * I->molar_weight[0] +
				concentration(I, 1, i, j, k) * I->molar_weight[1] +
				concentration(I, 2, i, j, k) * I->molar_weight[2] +
				concentration(I, 3, i, j, k) * I->molar_weight[3]) / (I->R * temperature_flow(I, i, j, k));
	else {
		printf("Error density derivative by pressure index\n");
		return 0;
	}
	nan_check(I, tmp, p, -1, i, j, k, __func__);
	return tmp;
}

double Darsi_A_coef_func(in *I, int i, int j, int k)
{
	double tmp = 0;
	for (int p = 0; p < 3; p++) {
		tmp += I->porousness[POR_IND(I, i, j, k)] * saturation(I, p, i, j, k) * density_derivative_by_pressure(I, p, i, j, k);
	}
	nan_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double Darsi_A_coef(in *I, int i, int j, int k)
{
	return I->array_of_parameters[PARAM_IND(I, 66, i, j, k)];
}

double thermal_conductivity_func(in *I, int i, int j, int k)
{
	double tmp = 0;
	for (int pp = 0; pp < 3; pp++)
		tmp += I->thermal_conductivity_coef[pp] * saturation(I, pp, i, j, k);
	nan_check(I, tmp, -1, -1, i, j, k, __func__);
	return tmp;
}

double thermal_conductivity(in *I, int i, int j, int k)
{
	return I->array_of_parameters[PARAM_IND(I, 67, i, j, k)];
}

double convection_speed_func(in *I, int pr, int i, int j, int  k)
{
	if (!((pr == 0) || (pr == 1) || (pr == 2)))
		printf("Error index of convection speed function\n");
	double tmp = 0, tmp1 = 0;
	for (int pp = 0; pp < 2; pp++)
		tmp += density_t(I, pp, i, j, k) * I->specific_heat[pp] * I->adiabatic_exponent[pp] * avarage_velocity(I, pp, pr, i, j, k) / (2 * I->dx[pr]);
	for (int pp = 0; pp < 4; pp++)
		tmp += density_t(I, 2, i, j, k) * I->specific_heat[pp + 2] * I->adiabatic_exponent[pp + 2] * concentration(I, pp, i, j, k) * avarage_velocity(I, 2, pr, i, j, k) / (2 * I->dx[pr]);
	for (int pp = 0; pp < 2; pp++)
		tmp1 += I->porousness[POR_IND(I, i, j, k)] * density_t(I, pp, i, j, k) * saturation(I, pp, i, j, k) * I->specific_heat[pp];
	for (int pp = 0; pp < 4; pp++)
		tmp1 += I->porousness[POR_IND(I, i, j, k)] * density_t(I, 2, i, j, k) * saturation(I, 2, i, j, k) * I->specific_heat[pp + 2] * concentration(I, pp, i, j, k);
	if (tmp1 == 0)
		return 0;
	else
		return tmp/tmp1;
}

double convection_speed(in *I, int pr, int i, int j, int  k)
{
	if (!((pr == 0) || (pr == 1) || (pr == 2)))
		printf("Error index of convection speed\n");
	return I->array_of_parameters[PARAM_IND(I, 68 + pr, i, j, k)];
}

int check_sum(in *I)
{
	int i, j, k, p, flag_too_small_saturation = 0, flag_too_big_saturation = 0, flag_bad_sum_saturation = 0;
	int flag_too_small_concentration = 0, flag_too_big_concentration = 0, flag_bad_sum_concentration = 0;
	int flag_too_small_temperature_flow = 0, flag_too_small_temperature_environment = 0;
	double residual_saturations_sum = 0, saturation_sum, concentration_sum;
	double disparity = 0;
	int disp_cells = 0;
	for (int p = 0; p < 3; p++)
		residual_saturations_sum += I->residual_saturation[p];
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
					disparity += fabs(1 - saturation_sum);
					if (fabs(1 - saturation_sum) > I->epsilon)
						disp_cells++;
					for (p = 0; p < 3; p++) {
						nan_check(I, saturation(I, p, i, j, k), p, -1, i, j, k, "saturation_before_normalization");
						if (saturation(I, p, i, j, k) < I->residual_saturation[p]) {
							I->B_prev[B_IND(I, p + 5, i, j, k)] = 0;
							flag_too_small_saturation = 1;
						}
						//if (saturation(I, p, i, j, k) > 1 - residual_saturations_sum + I->residual_saturation[p]) {
						//	I->B_prev[B_IND(I, p + 5, i, j, k)] = 1 - residual_saturations_sum;
						//	flag_too_big_saturation = 1;
						//}
					}
					saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
					if (saturation_sum > 1) {
						I->B_prev[B_IND(I, 5, i, j, k)] = (saturation(I, 0, i, j, k) - I->residual_saturation[0]) *
							(1 - saturation(I, 1, i, j, k) - I->residual_saturation[0] - I->residual_saturation[2]) /
							(saturation_sum - saturation(I, 1, i, j, k) - I->residual_saturation[0] - I->residual_saturation[2]);
						I->B_prev[B_IND(I, 7, i, j, k)] = (saturation(I, 2, i, j, k) - I->residual_saturation[2]) *
							(1 - saturation(I, 1, i, j, k) - I->residual_saturation[0] - I->residual_saturation[2]) /
							(saturation_sum - saturation(I, 1, i, j, k) - I->residual_saturation[0] - I->residual_saturation[2]);
						saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
						if (saturation(I, 0, i, j, k) < I->residual_saturation[0])
							I->B_prev[B_IND(I, 5, i, j, k)] = 0;
						if (saturation(I, 2, i, j, k) < I->residual_saturation[2])
							I->B_prev[B_IND(I, 7, i, j, k)] = 0;
						saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
						if (saturation_sum > 1)
							I->B_prev[B_IND(I, 6, i, j, k)] -= saturation_sum - 1;
						/*
						I->B_prev[B_IND(I, 6, i ,j, k)] -= saturation_sum - 1;
						if (saturation(I, 1, i, j, k) < I->residual_saturation[1]) {
							I->B_prev[B_IND(I, 6, i, j, k)] = 0;
							saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
							I->B_prev[B_IND(I, 7, i, j, k)] -= saturation_sum - 1;
							if (saturation(I, 2, i, j, k) < I->residual_saturation[2]) {
								I->B_prev[B_IND(I, 7, i, j, k)] = 0;
								saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
								I->B_prev[B_IND(I, 5, i, j, k)] -= saturation_sum - 1;
							}
						}
						*/
						flag_bad_sum_saturation = 1;
					} else if (saturation_sum < 1) {
						I->B_prev[B_IND(I, 6, i, j, k)] += 1 - saturation_sum;
						flag_bad_sum_saturation = 1;
					}
					saturation_sum = saturation(I, 0, i, j, k) + saturation(I, 1, i, j, k) + saturation(I, 2, i, j, k);
					if (fabs(saturation_sum - 1.0) > I->epsilon)
						printf("Error saturation sum after normalisation! i = %d, j = %d, k = %d, saturation_sum = %e\n", i, j, k, saturation_sum);
					for (p = 0; p < 3; p++)
						nan_check(I, saturation(I, p, i, j, k), p, -1, i, j, k, "saturation_after_normalization");
					for (p = 0, concentration_sum = 0; p < 4; p++) {
						nan_check(I, concentration(I, p, i, j, k), p, -1, i, j, k, "concentration_before_normalization");
						if (concentration(I, p, i, j, k) < 0) {
							I->B_prev[B_IND(I, p, i, j, k)] = 0;
							flag_too_small_concentration = 1;
						}
						//if (concentration(I, p, i, j, k) > 1) {
						//	I->B_prev[B_IND(I, p, i, j, k)] = 1;
						//	flag_too_big_concentration = 1;
						//}
						concentration_sum += concentration(I, p, i, j, k);
					}
					if (concentration_sum < I->epsilon) {
						I->B_prev[B_IND(I, 3, i, j, k)] = I->epsilon;
						concentration_sum = concentration(I, 0, i, j, k) + concentration(I, 1, i, j, k) + concentration(I, 2, i, j, k) + concentration(I, 3, i, j, k);
					}
					if (concentration_sum != 1) {
						I->B_prev[B_IND(I, 0, i, j, k)] = concentration(I, 0, i, j, k) / concentration_sum;
						I->B_prev[B_IND(I, 1, i, j, k)] = concentration(I, 1, i, j, k) / concentration_sum;
						I->B_prev[B_IND(I, 2, i, j, k)] = concentration(I, 2, i, j, k) / concentration_sum;
						I->B_prev[B_IND(I, 3, i, j, k)] = concentration(I, 3, i, j, k) / concentration_sum;
						flag_bad_sum_concentration = 1;
					}
					concentration_sum = concentration(I, 0, i, j, k) + concentration(I, 1, i, j, k) + concentration(I, 2, i, j, k) + concentration(I, 3, i, j, k);
					if (fabs(concentration_sum - 1.0) > I->epsilon)
						printf("Error concentration sum after normalisation! i = %d, j = %d, k = %d, concentration_sum = %e\n", i, j, k, concentration_sum);
					for (p = 0; p < 4; p++)
						nan_check(I, concentration(I, p, i, j, k), p, -1, i, j, k, "concentration_after_normalization");
					if (temperature_flow(I, i, j, k) < 400) {
						I->B_prev[B_IND(I, 8, i, j, k)] = 0;
						flag_too_small_temperature_flow = 1;
					}
					if (temperature_environment(I, i, j, k) < 400) {
						I->B_prev[B_IND(I, 9, i, j, k)] = 0;
						flag_too_small_temperature_environment = 1;
					}
				}
			}
		}
	}
	if (disparity != 0)
		printf("Disparity is %lf in %d cells\n", disparity, disp_cells);
	if (flag_too_small_saturation)
		printf("Check sum: There is too small saturation at time step %d\n", I->time_step);
	if (flag_too_big_saturation)
		printf("Check sum: There is too big saturation at time step %d\n", I->time_step);
	if (flag_too_small_concentration)
		printf("Check sum: There is too small concentration at time step %d\n", I->time_step);
	if (flag_too_big_concentration)
		printf("Check sum: There is too big concentration at time step %d\n", I->time_step);
	if (flag_bad_sum_saturation)
		printf("Check sum: There is bad saturations sum at time step %d\n", I->time_step);
	if (flag_bad_sum_concentration)
		printf("Check sum: There is bad concentrations sum at time step %d\n", I->time_step);
	if (flag_too_small_temperature_flow)
		printf("Check sum: There is too small temperature flow at time step %d\n", I->time_step);
	if (flag_too_small_temperature_environment)
		printf("Check sum: There is too small temperature environment at time step %d\n", I->time_step);
	return 0;
}

int check_convection_speed(in *I)
{
	double max_conv_speed_magnitude = 0, conv_speed_magnitude, max_grad_pressure_magnitude = 0, grad_pressure_magnitude;
	for (int k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (int i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (int j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if (cell_of_computation_domain(I, i, j, k)) {
					conv_speed_magnitude = sqrt(pow(convection_speed(I, 0, i, j, k) , 2) + pow(convection_speed(I, 1, i, j, k), 2));
					if (conv_speed_magnitude > max_conv_speed_magnitude)
						max_conv_speed_magnitude = conv_speed_magnitude;
					grad_pressure_magnitude = sqrt(pow(grad_pressure(I, 0, i, j, k) , 2) + pow(grad_pressure(I, 1, i, j, k), 2));
					if (grad_pressure_magnitude > max_grad_pressure_magnitude)
						max_grad_pressure_magnitude = grad_pressure_magnitude;
				}
			}
		}
	}
	if (max_conv_speed_magnitude > I->avarage_velocity_global_value) {
		printf("Too big convection speed!\n");
		max_conv_speed_magnitude = I->avarage_velocity_global_value;
	}
	for (int k = -I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (int i = -I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (int j = -I->stencil_size; j < I->ny + I->stencil_size; j++) {
				if (cell_of_computation_domain(I, i, j, k)) {
					I->array_of_parameters[PARAM_IND(I, 68, i, j, k)] = grad_pressure(I, 0, i, j, k) * max_conv_speed_magnitude / max_grad_pressure_magnitude;
					I->array_of_parameters[PARAM_IND(I, 69, i, j, k)] = grad_pressure(I, 1, i, j, k) * max_conv_speed_magnitude / max_grad_pressure_magnitude;
					I->array_of_parameters[PARAM_IND(I, 70, i, j, k)] = grad_pressure(I, 2, i, j, k) * max_conv_speed_magnitude / max_grad_pressure_magnitude;
				}
			}
		}
	}
	return 0;
}

int check_values_in_B(in *I)
{
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx; i++) {
			for (int j = 0; j < I->ny; j++) {
				for (int p = 0; p < I->num_parameters; p++) {
					if (I->B[A_IND(I, p, i, j, k)] != I->B[A_IND(I, p, i, j, k)]) {
						printf("Value %lf is written in B in (%d, %d, %d, %d)\n", I->B[A_IND(I, p, i, j, k)], p, i, j, k);
						I->nan_flag = 1;
					}
				}
			}
		}
	}
	return 0;
}

int check_values_in_B_prev(in *I)
{
	for (int k = - I->stencil_size; k < I->nz + I->stencil_size; k++) {
		for (int i = - I->stencil_size; i < I->nx + I->stencil_size; i++) {
			for (int j = - I->stencil_size; j < I->ny + I->stencil_size; j++) {
				for (int p = 0; p < I->num_parameters; p++) {
					if (I->B_prev[B_IND(I, p, i, j, k)] != I->B_prev[B_IND(I, p, i, j, k)]) {
						printf("Value %lf is written in B in (%d, %d, %d, %d)\n", I->B_prev[B_IND(I, p, i, j, k)], p, i, j, k);
						I->nan_flag = 1;
					}
				}
			}
		}
	}
	return 0;
}

int check_for_symmetry(in *I)
{
	for (int k = 0; k < I->nz; k++) {
		for (int i = 0; i < I->nx / 2; i++) {
			for (int j = 0; j <= i; j++) {
				for (int p = 0; p < I->num_parameters; p++) {
					if ((fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, I->nx - i, j, k)]) < 1e-10) ||
						(fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, i, I->ny - j, k)]) < 1e-10) ||
						(fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, I->nx - i, I->ny - j, k)]) < 1e-10) ||
						(fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, j, i, k)]) < 1e-10) ||
						(fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, I->ny - j, i, k)]) < 1e-10) ||
						(fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, j, I->nx - i, k)]) < 1e-10) ||
						(fabs(I->B_prev[B_IND(I, p, i, j, k)] - I->B_prev[B_IND(I, p, I->ny - j, I->nx - i, k)]) < 1e-10))
							printf("There is asymmetry in parameter %d in (%d,%d, %d)\n", p, i, j, k);
				}
			}
		}
	}
}

int print_oil_production(in *I)
{
	FILE *f;
	int production_well_coordinates[3];
	production_well_coordinates[0] = production_well_coordinates[1] = production_well_coordinates[2] = 0;
	int injection_well_coordinates[3];
	injection_well_coordinates[0] = I->gl_nx / 2;
	injection_well_coordinates[1] = I->gl_ny / 2;
	injection_well_coordinates[2] = 0;

	double velocity_of_oil_production =
		saturation(I, 1, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2]) *
		I->porousness[POR_IND(I, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2])] * I->dx[0] * I->dx[1] * I->dx[2] -
		I->B[A_IND(I, 6, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2])] *
		I->porousness[POR_IND(I, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2])] * I->dx[0] * I->dx[1] * I->dx[2] -
		avarage_velocity(I, 1, 0, production_well_coordinates[0] + 1, production_well_coordinates[1], production_well_coordinates[2]) * I->dx[1] * I->dx[2] * I->dt -
		avarage_velocity(I, 1, 1, production_well_coordinates[0], production_well_coordinates[1] + 1, production_well_coordinates[2]) * I->dx[0] * I->dx[2] * I->dt;
	if (I->time_step == 0) {
		if ((f = fopen("result/velocity_of_oil_production_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf",\
			I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness[POR_IND(I, 1, 1, 1)] * density_t(I, 1, 1, 1, 1));
	} else {
		if ((f = fopen("result/velocity_of_oil_production_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, velocity_of_oil_production / I->dt);
	fclose(f);

	I->volume_producted_oil_m += velocity_of_oil_production;
	if (I->time_step == 0) {
		if ((f = fopen("result/oil_production_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf\n%lf\t%lf\n",\
			I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness[POR_IND(I, 1, 1, 1)] * density_t(I, 1, 1, 1, 1), 0.0, 0.0);
	} else {
		if ((f = fopen("result/oil_production_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, I->volume_producted_oil_m);
	fclose(f);

	velocity_of_oil_production *= density_t(I, 1, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2]);
	if (I->time_step == 0) {
		if ((f = fopen("result/velocity_of_oil_production_kg.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf\n%lf\t%lf\n",\
			I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness[POR_IND(I, 1, 1, 1)] * density_t(I, 1, 1, 1, 1), 0.0, 0.0);
	} else {
		if ((f = fopen("result/velocity_of_oil_production_kg.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, velocity_of_oil_production / I->dt);
	fclose(f);

	I->volume_producted_oil_kg += velocity_of_oil_production;
	if (I->time_step == 0) {
		if ((f = fopen("result/oil_production_kg.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total oil volume %lf\n%lf\t%lf\n",\
			I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * saturation(I, 1, 1, 1, 1) * I->porousness[POR_IND(I, 1, 1, 1)] * density_t(I, 1, 1, 1, 1), 0.0, 0.0);
	} else {
		if ((f = fopen("result/oil_production_kg.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, I->volume_producted_oil_kg);
	fclose(f);

	double velocity_of_fluid_production = 0;
	for (int p = 0; p < 3; p++) {
		velocity_of_fluid_production +=
			saturation(I, p, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2]) *
			I->porousness[POR_IND(I, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2])] * I->dx[0] * I->dx[1] * I->dx[2] -
			I->B[A_IND(I, p + 5, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2])] *
			I->porousness[POR_IND(I, production_well_coordinates[0], production_well_coordinates[1], production_well_coordinates[2])] * I->dx[0] * I->dx[1] * I->dx[2] -
			avarage_velocity(I, p, 0, production_well_coordinates[0] + 1, production_well_coordinates[1], production_well_coordinates[2]) * I->dx[1] * I->dx[2] * I->dt -
			avarage_velocity(I, p, 1, production_well_coordinates[0], production_well_coordinates[1] + 1, production_well_coordinates[2]) * I->dx[0] * I->dx[2] * I->dt;
	}
	if (I->time_step == 0) {
		if ((f = fopen("result/velocity_of_fluid_production_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total fluid volume %lf\n%lf\t%lf\n",\
			I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * I->porousness[POR_IND(I, 1, 1, 1)], 0.0, 0.0);
	} else {
		if ((f = fopen("result/velocity_of_fluid_production_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, velocity_of_fluid_production / I->dt);
	fclose(f);

	I->volume_producted_fluid_m += velocity_of_fluid_production;
	if (I->time_step == 0) {
		if ((f = fopen("result/fluid_production_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
		fprintf(f, "#total fluid volume %lf\n%lf\t%lf\n", I->nx * I->dx[0] * I->ny * I->dx[1] * I->nz * I->dx[2] * I->porousness[POR_IND(I, 1, 1, 1)], 0.0, 0.0);
	} else {
		if ((f = fopen("result/fluid_production_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, I->volume_producted_fluid_m);
	fclose(f);

	double velocity_of_fluid_injection = 0;
	for (int p = 0; p < 3; p += 2) {
		velocity_of_fluid_injection +=
			avarage_velocity(I, p, 0, injection_well_coordinates[0] + 1, injection_well_coordinates[1], injection_well_coordinates[2]) * I->dx[1] * I->dx[2] * I->dt +
			avarage_velocity(I, p, 1, injection_well_coordinates[0], injection_well_coordinates[1] + 1, injection_well_coordinates[2]) * I->dx[0] * I->dx[2] * I->dt;
	}
	if (I->time_step == 0) {
		if ((f = fopen("result/velocity_of_fluid_injection_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
	} else {
		if ((f = fopen("result/velocity_of_fluid_injection_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, velocity_of_fluid_injection / I->dt);
	fclose(f);

	I->volume_injected_fluid_m += velocity_of_fluid_injection;
	if (I->time_step == 0) {
		if ((f = fopen("result/fluid_injection_m.dat","w")) == NULL) {
			printf("error openning file");
			return 1;
		}
	} else {
		if ((f = fopen("result/fluid_injection_m.dat","a")) == NULL) {
			printf("error openning file");
			return 1;
		}
	}
	fprintf(f, "%lf\t%lf\n", I->time, I->volume_injected_fluid_m);
	fclose(f);
	return 0;
}

double avarage_velocity_global(in *I)
{
	int i, j, k, p, pr;
	double y[4], tmp, min_saturation_gas = 1, max_abs_avarage_velocity_gas = 0;
	//printf("avaraje velocity global p = 0 pr = 0 i = 23 j = 24 %lf\n", avarage_velocity(I, 0, 0, 23, 24, 0) / (saturation(I, 2, 23, 24, 0) * I->porousness));
	//printf("avaraje velocity global p = 0 pr = 1 i = 23 j = 24 %lf\n", avarage_velocity(I, 0, 1, 23, 24, 0) / (saturation(I, 2, 23, 24, 0) * I->porousness));
	//printf("avaraje velocity global p = 1 pr = 0 i = 23 j = 24 %lf\n", avarage_velocity(I, 1, 0, 23, 24, 0) / (saturation(I, 2, 23, 24, 0) * I->porousness));
	//printf("avaraje velocity global p = 1 pr = 1 i = 23 j = 24 %lf\n", avarage_velocity(I, 1, 1, 23, 24, 0) / (saturation(I, 2, 23, 24, 0) * I->porousness));
	//printf("avaraje velocity global p = 2 pr = 0 i = 23 j = 24 %lf\n", avarage_velocity(I, 2, 0, 23, 24, 0) / (saturation(I, 2, 23, 24, 0) * I->porousness));
	//printf("avaraje velocity global p = 2 pr = 1 i = 23 j = 24 %lf\n", avarage_velocity(I, 2, 1, 23, 24, 0) / (saturation(I, 2, 23, 24, 0) * I->porousness));
	y[0] = y[1] = y[2] = y[3] = 0;
	for (k = 0; k < I->nz; k++) {
		for (i = 0; i < I->nx; i++) {
			for (j = 0; j < I->ny; j++) {
				if (I->ind_cell_multipl[i * I->ny + j] != -1) {
					for (p = 0; p < 3; p++) {
						tmp = 0;
						tmp = sqrt(pow(avarage_velocity(I, p, 0, i, j, k) / (saturation(I, p, i, j, k) * I->porousness[POR_IND(I, i, j, k)]), 2) +
							pow(avarage_velocity(I, p, 1, i ,j ,k) / (saturation(I, p, i, j, k) * I->porousness[POR_IND(I, i, j, k)]), 2));
						//for (pr = 0; pr < 2; pr++) {
						//	tmp += fabs(avarage_velocity(I, p, pr, i, j, k) / (saturation(I, p, i, j, k) * I->porousness[POR_IND(I, i, j, k)]));
						//}
						if (y[p] < tmp)
							y[p] = tmp;
					}
					if (min_saturation_gas > saturation(I, 2, i, j, k))
						min_saturation_gas = saturation(I, 2, i, j, k);
					if (max_abs_avarage_velocity_gas < fabs(avarage_velocity(I, 2, 0, i, j, k)) + fabs(avarage_velocity(I, 2, 1, i, j, k)))
						max_abs_avarage_velocity_gas = fabs(avarage_velocity(I, 2, 0, i, j, k)) + fabs(avarage_velocity(I, 2, 1, i, j, k));
				}
			}
		}
	}
	for (p = 0; p < 3; p++) {
		printf("Max_avarage_velocity%d = %lf\n", p, y[p]);
		if (y[3] < y[p])
			y[3] = y[p];
	}
	return y[3];
	if (min_saturation_gas == 0)
		min_saturation_gas = I->epsilon;
	if (max_abs_avarage_velocity_gas / (min_saturation_gas * I->porousness[POR_IND(I, i, j, k)]) > y[3])
		y[3] = max_abs_avarage_velocity_gas / (min_saturation_gas * I->porousness[POR_IND(I, i, j, k)]);
//	x /= I->n_cells_multipl * 9;
	return y[3];
}

int save_I_termogas(in *I)
{
	FILE *f = fopen("result/I", "w");
	fprintf(f, "%s\n", I->map_name);
	//fprintf(f, "%s\n", I->region_map_name);
	fprintf(f, "%lf\n", I->hight);
	fprintf(f, "%d %d %d\n", I->kx, I->ky, I->kz);
	fprintf(f, "%d %d %d\n", I->nx, I->ny, I->nz);
	fprintf(f, "%d\n", I->ncols);
	fprintf(f, "%d\n", I->nrows);
	fprintf(f, "%lf\n", I->cellsize);
	fprintf(f, "%d %d\n", I->n_points, I->n_cells);
	fprintf(f, "%d\n", I->n_points_multipl);
	fprintf(f, "%d\n", I->n_cells_multipl);
	//fprintf(f, "%lf %lf\n", I->interpolation, I->interpolation_poli);
	fprintf(f, "%lf\n", I->nodata_value);
	fprintf(f, "%lf %lf %lf\n", I->dx[0], I->dx[1], I->dx[2]);
	fprintf(f, "%d\n", I->num_parameters);
	fprintf(f, "%d %d\n", I->system_dimension, I->system_dimension_with_boundary);
	fprintf(f, "%d\n", I->non_zero_elem);
	fprintf(f, "%lf\n", I->end_time);
	fprintf(f, "%d\n", I->stencil_size);
	fprintf(f, "%d\n", I->n_boundary_cells);
	fprintf(f, "%lf\n", I->mass_quantity);
	//fprintf(f, "%lf\n", I->porousness);
	fprintf(f, "%lf %lf\n", I->density_0[0], I->density_0[1]);
	fprintf(f, "%lf\n", I->pressure_0);
	fprintf(f, "%lf %lf\n", I->density_coef_beta[0], I->density_coef_beta[1]);
	fprintf(f, "%lf\n", I->temperature_0);
	fprintf(f, "%lf %lf\n", I->density_coef_a[0], I->density_coef_a[1]);
	fprintf(f, "%lf\n", I->R);
	fprintf(f, "%lf %lf %lf %lf %lf\n", I->molar_weight[0],
		I->molar_weight[1], I->molar_weight[2], I->molar_weight[3],
		I->molar_weight[4]);
	fprintf(f, "%lf\n", I->permeability);
	fprintf(f, "%lf %lf %lf\n", I->residual_saturation[0], I->residual_saturation[1], I->residual_saturation[2]);
	fprintf(f, "%lf %lf\n", I->viscosity_coef_A[0], I->viscosity_coef_A[1]);
	fprintf(f, "%lf %lf %lf %lf\n", I->viscosity_coef_A_gas[0], I->viscosity_coef_A_gas[1],
		I->viscosity_coef_A_gas[2], I->viscosity_coef_A_gas[3]);
	fprintf(f, "%lf %lf\n", I->viscosity_coef_B[0], I->viscosity_coef_B[1]);
	fprintf(f, "%lf %lf %lf %lf\n", I->viscosity_coef_C_gas[0], I->viscosity_coef_C_gas[1],
		I->viscosity_coef_C_gas[2], I->viscosity_coef_C_gas[3]);
	fprintf(f, "%lf %lf %lf %lf\n", I->temperature_0_gas[0], I->temperature_0_gas[1],
		I->temperature_0_gas[2], I->temperature_0_gas[3]);
	fprintf(f, "%lf %lf %lf\n", I->capillary_pressure_at_maximum_saturation[0],
		I->capillary_pressure_at_maximum_saturation[1], I->capillary_pressure_at_maximum_saturation[2]);
	fprintf(f, "%lf %lf %lf\n", I->residual_saturation_two_phase[0],
		I->residual_saturation_two_phase[1], I->residual_saturation_two_phase[2]);
	fprintf(f, "%lf\n", I->capillary_pressure_coef);
	fprintf(f, "%lf\n", I->num_carbon_atoms);
	fprintf(f, "%lf\n", I->num_hydrogen_atoms);
	fprintf(f, "%lf %lf %lf\n", I->stoichiometric_coef[0], I->stoichiometric_coef[1], I->stoichiometric_coef[2]);
	fprintf(f, "%lf %lf %lf %lf %lf %lf\n", I->stoichiometric_coef_before[0], I->stoichiometric_coef_before[1],
		I->stoichiometric_coef_before[2], I->stoichiometric_coef_before[3], I->stoichiometric_coef_before[4],
		I->stoichiometric_coef_before[5]);
	fprintf(f, "%lf %lf %lf %lf %lf %lf\n", I->stoichiometric_coef_after[0], I->stoichiometric_coef_after[1],
		I->stoichiometric_coef_after[2], I->stoichiometric_coef_after[3], I->stoichiometric_coef_after[4],
		I->stoichiometric_coef_after[5]);
	fprintf(f, "%lf\n", I->threshold_temperature);
	fprintf(f, "%lf\n", I->frequency_factor);
	fprintf(f, "%lf\n", I->activation_temperature);
	fprintf(f, "%lf\n", I->stoichiometric_coef_activ);
	fprintf(f, "%lf\n", I->pressure_activ);
	fprintf(f, "%lf %lf %lf %lf %lf %lf %lf\n", I->specific_heat[0], I->specific_heat[1],
		I->specific_heat[2], I->specific_heat[3], I->specific_heat[4], I->specific_heat[5], I->specific_heat[6]);
	fprintf(f, "%lf %lf %lf %lf\n", I->thermal_conductivity_coef[0], I->thermal_conductivity_coef[1],
		I->thermal_conductivity_coef[2], I->thermal_conductivity_coef[3]);
	fprintf(f, "%lf\n", I->heat_transfer_coef);
	fprintf(f, "%lf\n", I->tempetarure_for_calculation_internal_energy);
	fprintf(f, "%lf %lf %lf %lf %lf %lf %lf\n", I->initial_enthalpy[0], I->initial_enthalpy[1],
		I->initial_enthalpy[2], I->initial_enthalpy[3], I->initial_enthalpy[4], I->initial_enthalpy[5], I->initial_enthalpy[6]);
	fprintf(f, "%lf\n", I->density_environment);
	fprintf(f, "%lf\n", I->initial_pressure);
	fprintf(f, "%lf\n", I->initial_temperature);
	fprintf(f, "%lf\n", I->injection_well_pressure);
	fprintf(f, "%lf\n", I->injection_well_temperature);
	fprintf(f, "%lf\n", I->production_well_pressure);
	fprintf(f, "%d\n", I->equation_num);
	fprintf(f, "%lf\n", I->epsilon);
	fprintf(f, "%lf %lf %lf %lf %lf %lf\n", I->adiabatic_exponent[0], I->adiabatic_exponent[1],
		I->adiabatic_exponent[2], I->adiabatic_exponent[3], I->adiabatic_exponent[4], I->adiabatic_exponent[5]);
	fprintf(f, "%d\n", I->time_step);
	fprintf(f, "%lf\n", I->volume_producted_oil_kg);
	fprintf(f, "%lf\n", I->volume_producted_oil_m);
	fprintf(f, "%lf\n", I->volume_producted_fluid_m);
	fprintf(f, "%lf\n", I->volume_injected_fluid_m);
	fprintf(f, "%d\n", I->my_rank);
	fprintf(f, "%d\n", I->nproc);
	fprintf(f, "%d %d %d\n", I->gl_nx, I->gl_ny, I->gl_nz);
	fprintf(f, "%d %d\n", I->x_regions, I->y_regions);
	fprintf(f, "%d %d\n", I->num_el_in_x_region, I->num_el_in_y_region);
	fprintf(f, "%d %d\n", I->max_num_el_in_x_region, I->max_num_el_in_y_region);
	fprintf(f, "%d\n", I->gl_n_cells_multipl);
	fclose(f);
	return 0;
}

int read_I_termogas(in *I)
{
	FILE *f = fopen("result/I", "r");
	char s[100];
	fgets(s, 100, f);
	sscanf(s, "%s\n", I->map_name);
	fgets(s, 100, f);
	sscanf(s, "%s\n", I->region_map_name);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->hight);
	fgets(s, 100, f);
	sscanf(s, "%d %d %d\n", &I->kx, &I->ky, &I->kz);
	fgets(s, 100, f);
	sscanf(s, "%d %d %d\n", &I->nx, &I->ny, &I->nz);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->ncols);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->nrows);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->cellsize);
	fgets(s, 100, f);
	sscanf(s, "%d %d\n", &I->n_points, &I->n_cells);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->n_points_multipl);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->n_cells_multipl);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf\n", &I->interpolation, &I->interpolation_poli);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->nodata_value);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf\n", &I->dx[0], &I->dx[1], &I->dx[2]);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->num_parameters);
	fgets(s, 100, f);
	sscanf(s, "%d %d\n", &I->system_dimension, &I->system_dimension_with_boundary);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->non_zero_elem);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->end_time);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->stencil_size);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->n_boundary_cells);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->mass_quantity);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->porousness);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf\n", &I->density_0[0], &I->density_0[1]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->pressure_0);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf\n", &I->density_coef_beta[0], &I->density_coef_beta[1]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->temperature_0);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf\n", &I->density_coef_a[0], &I->density_coef_a[1]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->R);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf %lf\n", &I->molar_weight[0],
				&I->molar_weight[1], &I->molar_weight[2], &I->molar_weight[3],
				&I->molar_weight[4]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->permeability);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf\n", &I->residual_saturation[0], &I->residual_saturation[1], &I->residual_saturation[2]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf\n", &I->viscosity_coef_A[0], &I->viscosity_coef_A[1]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf\n", &I->viscosity_coef_A_gas[0], &I->viscosity_coef_A_gas[1],
	 	  	&I->viscosity_coef_A_gas[2], &I->viscosity_coef_A_gas[3]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf\n", &I->viscosity_coef_B[0], &I->viscosity_coef_B[1]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf\n", &I->viscosity_coef_C_gas[0], &I->viscosity_coef_C_gas[1],
		   	&I->viscosity_coef_C_gas[2], &I->viscosity_coef_C_gas[3]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf\n", &I->temperature_0_gas[0], &I->temperature_0_gas[1],
		   	&I->temperature_0_gas[2], &I->temperature_0_gas[3]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf\n", &I->capillary_pressure_at_maximum_saturation[0],
		   	&I->capillary_pressure_at_maximum_saturation[1], &I->capillary_pressure_at_maximum_saturation[2]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf\n", &I->residual_saturation_two_phase[0],
	 	  	&I->residual_saturation_two_phase[1], &I->residual_saturation_two_phase[2]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->capillary_pressure_coef);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->num_carbon_atoms);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->num_hydrogen_atoms);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf\n", &I->stoichiometric_coef[0], &I->stoichiometric_coef[1], &I->stoichiometric_coef[2]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf %lf %lf\n", &I->stoichiometric_coef_before[0], &I->stoichiometric_coef_before[1],
		   	&I->stoichiometric_coef_before[2], &I->stoichiometric_coef_before[3], &I->stoichiometric_coef_before[4],
		   	&I->stoichiometric_coef_before[5]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf %lf %lf\n", &I->stoichiometric_coef_after[0], &I->stoichiometric_coef_after[1],
		   	&I->stoichiometric_coef_after[2], &I->stoichiometric_coef_after[3], &I->stoichiometric_coef_after[4],
		   	&I->stoichiometric_coef_after[5]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->threshold_temperature);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->frequency_factor);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->activation_temperature);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->stoichiometric_coef_activ);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->pressure_activ);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf %lf %lf %lf\n", &I->specific_heat[0], &I->specific_heat[1],
		   	&I->specific_heat[2], &I->specific_heat[3], &I->specific_heat[4], &I->specific_heat[5], &I->specific_heat[6]);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf\n", &I->thermal_conductivity_coef[0], &I->thermal_conductivity_coef[1],
		   	&I->thermal_conductivity_coef[2], &I->thermal_conductivity_coef[3]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->heat_transfer_coef);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->tempetarure_for_calculation_internal_energy);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf %lf %lf %lf\n", &I->initial_enthalpy[0], &I->initial_enthalpy[1],
		   	&I->initial_enthalpy[2], &I->initial_enthalpy[3], &I->initial_enthalpy[4], &I->initial_enthalpy[5], &I->initial_enthalpy[6]);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->density_environment);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->initial_pressure);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->initial_temperature);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->injection_well_pressure);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->injection_well_temperature);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->production_well_pressure);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->equation_num);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->epsilon);
	fgets(s, 100, f);
	sscanf(s, "%lf %lf %lf %lf %lf %lf\n", &I->adiabatic_exponent[0], &I->adiabatic_exponent[1],
		   	&I->adiabatic_exponent[2], &I->adiabatic_exponent[3], &I->adiabatic_exponent[4], &I->adiabatic_exponent[5]);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->time_step);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->volume_producted_oil_kg);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->volume_producted_oil_m);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->volume_producted_fluid_m);
	fgets(s, 100, f);
	sscanf(s, "%lf\n", &I->volume_injected_fluid_m);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->my_rank);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->nproc);
	fgets(s, 100, f);
	sscanf(s, "%d %d %d\n", &I->gl_nx, &I->gl_ny, &I->gl_nz);
	fgets(s, 100, f);
	sscanf(s, "%d %d\n", &I->x_regions, &I->y_regions);
	fgets(s, 100, f);
	sscanf(s, "%d %d\n", &I->num_el_in_x_region, &I->num_el_in_y_region);
	fgets(s, 100, f);
	sscanf(s, "%d %d\n", &I->max_num_el_in_x_region, &I->max_num_el_in_y_region);
	fgets(s, 100, f);
	sscanf(s, "%d\n", &I->gl_n_cells_multipl);
	fclose(f);
	return 0;
}

int calculate_disparity_pressure(in *I)
{
	int x = 0;
	for(int i = 0; i < I->nx; i++) {
		for (int j = 0; j < I->ny; j++) {
			if (!well(I, i, j, 0)) {
				double tmp = 0.5 * ((Darsi_M_coef(I, i + 1, j, 0) + Darsi_M_coef(I, i, j, 0)) * (pressure(I, i + 1, j, 0) - pressure(I, i, j, 0)) -
					(Darsi_M_coef(I, i, j, 0) + Darsi_M_coef(I, i - 1, j, 0) * (pressure(I, i, j, 0) - pressure(I, i - 1, j, 0)))) / (I->dx[0] * I->dx[0]) +
					0.5 * ((Darsi_M_coef(I, i, j + 1, 0) + Darsi_M_coef(I, i, j, 0)) * (pressure(I, i, j + 1, 0) - pressure(I, i, j, 0)) -
					(Darsi_M_coef(I, i, j, 0) + Darsi_M_coef(I, i, j - 1, 0)) * (pressure(I, i, j, 0) - pressure(I, i, j - 1, 0))) / (I->dx[1] * I->dx[1]);
				if (tmp > 1e-03) {
					printf("Disparity of pressure in [%d, %d] = %e\n", i, j, tmp);
					x = 1;
				}
			}
		}
	}
	return x;
}

#endif
