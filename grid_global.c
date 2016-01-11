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
#include "grid.h"

#include "grid_global.h"


grid_t *initGrid_withDomain(domain_t *thisDomain, float boxsize, int gridsize)
{
	grid_t *newGrid;
	newGrid = initGrid();
	
	newGrid->gridsize = gridsize;
	
	for(int i=0; i<3; i++)
	{
		newGrid->upLimit[i] = thisDomain->uplimit[i]*boxsize;
		newGrid->lowLimit[i] = thisDomain->lowlimit[i]*boxsize;
		newGrid->upLimit_int[i] = roundf(thisDomain->uplimit[i]*gridsize);
		newGrid->lowLimit_int[i] = roundf(thisDomain->lowlimit[i]*gridsize);
		newGrid->dim[i] = (int)roundf(gridsize/thisDomain->dims[i]);
	}
	
	newGrid->clump = allocateGrid_withDomain_float(thisDomain, gridsize);
	newGrid->rho = allocateGrid_withDomain_float(thisDomain, gridsize);
	newGrid->inv_rho = allocateGrid_withDomain_float(thisDomain, gridsize);
	newGrid->npart_cell = allocateGrid_withDomain_int(thisDomain, gridsize);
	
	return newGrid;
}

float *allocateGrid_withDomain_float(domain_t *thisDomain, int gridsize)
{
	float *thisArray;
	int count;
	int dim[3];
	
	count = 1;
	for(int i=0; i<3; i++)
	{
		dim[i] = (int)roundf(gridsize/thisDomain->dims[i]);
		count *= dim[i];
	}
	
	thisArray = malloc(count*sizeof(float));
	if(thisArray == NULL)
	{
		fprintf(stderr, "ERROR: allocating array: Not enough memory to allocate array.\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<count; i++)
	{
		thisArray[i] = 0.f;
	}
	
	return thisArray;
}

int *allocateGrid_withDomain_int(domain_t *thisDomain, int gridsize)
{
	int *thisArray;
	int count;
	int dim[3];
	
	count = 1;
	for(int i=0; i<3; i++)
	{
		dim[i] = (int)roundf(gridsize/thisDomain->dims[i]);
		count *= dim[i];
	}
	
	thisArray = malloc(count*sizeof(int));
	if(thisArray == NULL)
	{
		fprintf(stderr, "ERROR: allocating array: Not enough memory to allocate array.\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<count; i++)
	{
		thisArray[i] = 0;
	}
	
	return thisArray;
}

void save_rho_to_file(grid_t *thisGrid, input_t *thisInput)
{
	char filename[200];
	  
	sprintf(filename, "%s_rho.dat", thisInput->fname);
	
  	int count = 1;
	
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
#endif

#ifdef __MPI
	MPI_Datatype mpi_subarray;
	
	int sizes[3], subsizes[3], starts[3];
	
	for(int i=0; i<3; i++)
	{
		sizes[i] = thisGrid->gridsize;
		subsizes[i] = thisGrid->dim[i];
		starts[i] = thisGrid->lowLimit_int[i];
		count *= subsizes[i];
	}
	MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &mpi_subarray);
	MPI_Type_commit(&mpi_subarray);
	
	offset = 0;
	MPI_File_set_view(mpifile, offset, MPI_FLOAT, mpi_subarray, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpifile, thisGrid->rho, count, MPI_FLOAT, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	
	
// 	float *buf;
// 	float *recv_buf;
// 	int skip, skip2;
// 	int starts_rec[3];
// 	int subsizes_rec[3];
// 	
// 	if(thisDomain->originRank != 0)
// 	{
// 		MPI_Ssend(&starts, 3, MPI_INT, 0, 101, MPI_COMM_WORLD);
// 		MPI_Ssend(&subsizes, 3, MPI_INT, 0, 102, MPI_COMM_WORLD);
// 		MPI_Ssend(thisGrid->rho, count, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
// 	}else{
// 		for(int sendRank=1; sendRank<thisDomain->size; sendRank++)
// 		{
// 			MPI_Recv(&starts_rec, 3, MPI_INT, sendRank, 101, MPI_COMM_WORLD, &status);
// 			MPI_Recv(&subsizes_rec, 3, MPI_INT, sendRank, 102, MPI_COMM_WORLD, &status);
// 			count = subsizes_rec[0]*subsizes[1]*subsizes[2];
// 			recv_buf = malloc(count*sizeof(float));
// 			MPI_Recv(recv_buf, count, MPI_FLOAT, sendRank, 100, MPI_COMM_WORLD, &status);
// 			
// 			offset = starts_rec[0]*thisGrid->gridsize*thisGrid->gridsize + starts_rec[1]*thisGrid->gridsize + starts_rec[2];
// 			skip = (thisGrid->gridsize-subsizes_rec[2]);
// 			skip2 = (thisGrid->gridsize*(thisGrid->gridsize-subsizes_rec[1]));
// 			for(int i=0, j=0; i<count; i++, j++)
// 			{
// 				if(i%(subsizes_rec[2]*subsizes_rec[1]) == 0 && i!=0) j = j + skip2;
// 				if(i%subsizes_rec[2] == 0 && i!=0) j = j+skip;
// 				buf[offset+ j] = recv_buf[i]*1000*sendRank;
// 
// 			}
// 			free(recv_buf);
// 		}
// 		for(int i=0, j=0; i<count; i++, j++)
// 		{
// 		  	skip = (thisGrid->gridsize-subsizes[2]);
// 			skip2 = (thisGrid->gridsize*(thisGrid->gridsize-subsizes[1]));
// 			if(i%(subsizes[2]*subsizes[1]) == 0 && i!=0) j = j + skip2;
// 			if(i%subsizes[2] == 0 && i!=0) j = j+skip;
// 			buf[j] = 0.;//thisGrid->rho[i];
// 		}
// 	}
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	
// 	if(thisDomain->originRank == 0)
// 	{
// 		FILE * fp;
// 		
// 		fp = fopen(filename, "wb");
// 		count = thisGrid->gridsize*thisGrid->gridsize*thisGrid->gridsize;
// 		fwrite(buf, sizeof(float), thisGrid->gridsize*thisGrid->gridsize*thisGrid->gridsize, fp);
// 		fclose(fp);
// 		
// 		free(buf);
// 	}
#else
	count = thisGrid->gridsize*thisGrid->gridsize*thisGrid->gridsize;
	fwrite(thisGrid->rho, sizeof(float), count, fp);
#endif

#ifdef __MPI
	MPI_File_close(&mpifile);
#else
	fclose(fp);
#endif
}

void save_inv_rho_to_file(grid_t *thisGrid, input_t *thisInput)
{
	char filename[200];
	  
	sprintf(filename, "%s_inv_rho.dat", thisInput->fname);
	
  	int count = 1;
	
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
#endif

#ifdef __MPI
	MPI_Datatype mpi_subarray;
	
	int sizes[3], subsizes[3], starts[3];
	
	for(int i=0; i<3; i++)
	{
		sizes[i] = thisGrid->gridsize;
		subsizes[i] = thisGrid->dim[i];
		starts[i] = thisGrid->lowLimit_int[i];
		count *= subsizes[i];
	}
	MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &mpi_subarray);
	MPI_Type_commit(&mpi_subarray);
	
	offset = 0;
	
	MPI_File_set_view(mpifile, offset, MPI_FLOAT, mpi_subarray, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpifile, thisGrid->inv_rho, count, MPI_FLOAT, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	
#else
	count = thisGrid->gridsize*thisGrid->gridsize*thisGrid->gridsize;
	fwrite(thisGrid->inv_rho, sizeof(float), count, fp);
#endif

#ifdef __MPI
	MPI_File_close(&mpifile);
#else
	fclose(fp);
#endif
}

void save_npart_cell_to_file(grid_t *thisGrid, input_t *thisInput)
{
	char filename[200];
	  
	sprintf(filename, "%s_npart_cells.dat", thisInput->fname);
	
  	int count = 1;
	
#ifdef __MPI
	MPI_File mpifile;
	MPI_Offset offset;
	MPI_Status status;
	
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
#else
	FILE * fp;
	
	fp = fopen(filename, "wb");
#endif

#ifdef __MPI
	MPI_Datatype mpi_subarray;
	
	int sizes[3], subsizes[3], starts[3];
	
	for(int i=0; i<3; i++)
	{
		sizes[i] = thisGrid->gridsize;
		subsizes[i] = thisGrid->dim[i];
		starts[i] = thisGrid->lowLimit_int[i];
		count *= subsizes[i];
	}
	MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &mpi_subarray);
	MPI_Type_commit(&mpi_subarray);
	
	offset = 0;
	
	MPI_File_set_view(mpifile, offset, MPI_FLOAT, mpi_subarray, "native", MPI_INFO_NULL);
	MPI_File_write_all(mpifile, thisGrid->npart_cell, count, MPI_FLOAT, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	
#else
	count = thisGrid->gridsize*thisGrid->gridsize*thisGrid->gridsize;
	fwrite(thisGrid->npart_cell, sizeof(int), count, fp);
#endif

#ifdef __MPI
	MPI_File_close(&mpifile);
#else
	fclose(fp);
#endif
}
