#ifndef GRID_H
#define GRID_H
#endif

typedef struct
{
	int gridsize;
	float *clump;
	
	float *rho;
	float *inv_rho;
	int *npart_cell;
	
	float upLimit[3];
	float lowLimit[3];
	int upLimit_int[3];
	int lowLimit_int[3];
	int dim[3];
} grid_t;

grid_t *initGrid();
void deallocateGrid(grid_t *thisGrid);
int totNcells_grid(grid_t * thisGrid);
float meanRho_grid(grid_t *thisGrid);

#ifdef __MPI
int totNcells_global(grid_t *thisGrid);
float meanRho_global(grid_t *thisGrid);
#endif