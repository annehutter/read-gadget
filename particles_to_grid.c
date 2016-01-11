#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "input.h"
#include "header.h"
#include "domain_decomposition.h"
#include "part.h"
#include "grid.h"

#include "grid_global.h"
#include "particles_to_grid.h"

double meanRho_part(part_t *theseParticles)
{
	double meanRho_part_local = 0;
	
	for(int p=0; p<theseParticles->num; p++)
	{
		meanRho_part_local += theseParticles->rho[p];
	}
	
	meanRho_part_local /= theseParticles->num;
	
	return meanRho_part_local;
}

double meanRho_part_global(part_t *theseParticles)
{
	int totNpart_global;
	double meanRho_part_local, meanRho_part_global;
	
	meanRho_part_local = meanRho_part(theseParticles)*theseParticles->num;
	MPI_Allreduce(&meanRho_part_local, &meanRho_part_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&theseParticles->num, &totNpart_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	meanRho_part_global /= totNpart_global;
	
	return meanRho_part_global;
}

void comp_arrays_for_clump(domain_t *thisDomain, part_t *theseParticles, grid_t *thisGrid, float boxsize)
{
	float upLimit = 100.;
  
	float inv_BoxSize = 1.f/boxsize;
	int gridsize = thisGrid->gridsize;
	float factor = inv_BoxSize*gridsize;
	int totNcells;
	int x_int, y_int, z_int;
	int index;
	float rho;
	float meanRhoPart;
	
	totNcells = totNcells_grid(thisGrid);
	
	printf("rank %d: totNCells = %d\n", thisDomain->originRank, totNcells);
	
	/* compute mean density of a particle */
	meanRhoPart = meanRho_part_global(theseParticles);
	
	printf("rank %d: meanRho_part = %e\n", thisDomain->originRank, meanRhoPart);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* go through particles in each domain and compute values */
	printf("rank %d: factor = %e\n", thisDomain->originRank, factor);
	printf("rank %d: low_limit = %d %d %d\n", thisDomain->originRank, thisGrid->lowLimit_int[0], thisGrid->lowLimit_int[1], thisGrid->lowLimit_int[2]);
	for(int p = 0; p<theseParticles->num; p++)
	{
		x_int = theseParticles->pos[3*p]*factor-thisGrid->lowLimit_int[0];
		y_int = theseParticles->pos[3*p+1]*factor-thisGrid->lowLimit_int[1];
		z_int = theseParticles->pos[3*p+2]*factor-thisGrid->lowLimit_int[2];
// 		printf("rank %d: %e %e %e: %d %d %d\t %d %d %d\t %d %d %d\t %e %e %e\n", thisDomain->originRank, theseParticles->pos[3*p], theseParticles->pos[3*p+1], theseParticles->pos[3*p+2], x_int, y_int, z_int, thisGrid->lowLimit_int[0], thisGrid->lowLimit_int[1], thisGrid->lowLimit_int[2], thisGrid->upLimit_int[0], thisGrid->upLimit_int[1], thisGrid->upLimit_int[2], thisDomain->lowlimit[0], thisDomain->lowlimit[1], thisDomain->lowlimit[2]);
		
		index = x_int*(thisGrid->dim[1]*thisGrid->dim[2]) + y_int*thisGrid->dim[2] + z_int;	//TODO: adapt to domain sized grid
		
		rho = theseParticles->rho[p];
		
		if(index >=0 && index < totNcells)
		{
			if(rho>0. && rho<upLimit*meanRhoPart)
			{
				thisGrid->rho[index] += rho;
				thisGrid->inv_rho[index] += 1./rho;
				thisGrid->npart_cell[index] += 1;
			}
		}else{
			printf("rank %d: particle %d has an incorrect index %d: %d %d %d\n", thisDomain->originRank, p, index, x_int, y_int, z_int);
		}
	}
}

void produce_clumping_factor_fields(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, int gridsize)
{
	float boxsize;
	if(thisDomain->originRank == 0)
	{
		boxsize = thisHeader->BoxSize;
	}
	MPI_Bcast(&boxsize, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	printf("rank %d: boxsize = %e\t %e\n",thisDomain->originRank, boxsize, thisHeader->BoxSize);
	
	grid_t *thisGrid;
	thisGrid = initGrid_withDomain(thisDomain, boxsize, gridsize);

	printf("rank %d: initialized grid\n", thisDomain->originRank);
	comp_arrays_for_clump(thisDomain, theseParticles, thisGrid, boxsize);
	
	save_rho_to_file(thisDomain, thisGrid, thisInput);
// 	save_rho_to_file(thisGrid, thisInput);
// 	save_rho_to_file(thisGrid, thisInput);
	
	deallocateGrid(thisGrid);
}