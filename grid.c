#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "grid.h"

grid_t *initGrid()
{
	grid_t *newGrid;
	newGrid = malloc(sizeof(grid_t));
	if(newGrid == NULL)
	{
		fprintf(stderr, "ERROR: initGrid: Not enough memory to allocate grid.\n");
		exit(EXIT_FAILURE);
	}
	
	newGrid->gridsize = 0;
	newGrid->clump = NULL;
	newGrid->rho = NULL;
	newGrid->inv_rho = NULL;
	newGrid->npart_cell = NULL;
	
	for(int i=0; i<3; i++)
	{
		newGrid->upLimit[i] = 0.f;
		newGrid->lowLimit[i] = 0.f;
		newGrid->upLimit_int[i] = 0;
		newGrid->lowLimit_int[i] = 0;
		newGrid->dim[i] = 0;
	}
	
	return newGrid;
}

void deallocateGrid(grid_t *thisGrid)
{
	if(thisGrid != NULL)
	{
		if(thisGrid->clump != NULL) free(thisGrid->clump);
		if(thisGrid->rho != NULL) free(thisGrid->rho);
		if(thisGrid->inv_rho != NULL) free(thisGrid->inv_rho);
		if(thisGrid->npart_cell != NULL) free(thisGrid->npart_cell);
		
		free(thisGrid);
	}
}

int totNcells_grid(grid_t *thisGrid)
{
	int totNcells = 1;
	
	for(int i=0; i<3; i++)
	{
		totNcells *= thisGrid->dim[i]; 
	}
	
	return totNcells;
}

float meanRho_grid(grid_t *thisGrid)
{
	float meanRho = 0.f;
	int totNcells = totNcells_grid(thisGrid);
	
	for(int i=0; i<totNcells; i++)
	{
		meanRho += thisGrid->rho[i];
	}
	
	return meanRho;
}

#ifdef __MPI
int totNcells_global(grid_t *thisGrid)
{
	int totNcells_local, totNcells;
	
	totNcells_local = totNcells_grid(thisGrid);
	MPI_Allreduce(&totNcells_local, &totNcells, 1, MPI_INT, MPI_SUM,MPI_COMM_WORLD);
	
	return totNcells;
}
#endif

#ifdef __MPI
float meanRho_global(grid_t *thisGrid)
{
	int totNcells;
	float meanRho_local, meanRho;
	
	meanRho_local = meanRho_grid(thisGrid);
	MPI_Allreduce(&meanRho_local, &meanRho, 1, MPI_FLOAT, MPI_SUM,MPI_COMM_WORLD);
	
	totNcells = totNcells_global(thisGrid);
	
	meanRho /= totNcells;
	
	return meanRho;
}
#endif