 #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "domain_decomposition.h"

domain_t *initDomain(int thisRank, int size)
{
	int *n;
	float tmpDim;
	
#ifdef __MPI
// 	MPI_Status status;
// 	MPI_Request request[2];
	MPI_Comm comm_cart;
	int mpiStat;

	
	int neighbourRank;
	int neighbourCoords[3];
	int exist = 0;
	int counter;
#endif
	
	int ndims = 3;
  
	domain_t *newDomain;
	newDomain = malloc(sizeof(domain_t));
	
	newDomain->ndims = ndims;
	
	if(((size%2)%3 != 0) && ((size%3)%2 !=0) && (size != 1))
	{
		fprintf(stderr, "ERROR: This number of processors is not allowed for domain decomposition.\n");
		exit(EXIT_FAILURE);
	}

	for(int i=0; i<ndims; i++)
	{
		newDomain->dims[i] = 1;
		newDomain->periods[i] = 1;
		newDomain->coords[i] = 0;
	}
	
	n = getNumCuts(size);
	getDims(ndims, n, newDomain->dims);
	
	newDomain->originRank = thisRank;
	newDomain->size = size;
	
	newDomain->numNeighbours = getNumNeighbours(size, ndims, n);
	
	if(newDomain->numNeighbours>0)
	{
		newDomain->listNeighbourRanks = malloc(newDomain->numNeighbours*sizeof(int));
		for(int i=0; i<newDomain->numNeighbours; i++)
		{
			newDomain->listNeighbourRanks[i] = -99;	//set to value, so rank 0 is not excluded later from being a neighbour
		}
		newDomain->neighbourCoords = malloc(newDomain->numNeighbours*sizeof(int)*3);
	}else{
		newDomain->listNeighbourRanks = NULL;
		newDomain->neighbourCoords = NULL;
	}
	
	if(size > 1)
	{
#ifdef __MPI
		mpiStat = MPI_Cart_create(MPI_COMM_WORLD, ndims, newDomain->dims, newDomain->periods, 1, &comm_cart);
		assert(mpiStat == MPI_SUCCESS);
		mpiStat = MPI_Cart_coords(comm_cart, thisRank, ndims, newDomain->coords);
		assert(mpiStat == MPI_SUCCESS);
	
		counter=0;
		for(int x=newDomain->coords[0]-1; x<=newDomain->coords[0]+1; x++){
			for(int y=newDomain->coords[1]-1; y<=newDomain->coords[1]+1; y++){
				for(int z=newDomain->coords[2]-1; z<=newDomain->coords[2]+1; z++){
					neighbourCoords[0]=x;
					neighbourCoords[1]=y;
					neighbourCoords[2]=z;
					mpiStat = MPI_Cart_rank(comm_cart,neighbourCoords,&neighbourRank);
					assert(mpiStat == MPI_SUCCESS);
					if(neighbourRank!=thisRank){
						exist = 0;
						for(int i=0; i<newDomain->numNeighbours; i++){
							if(newDomain->listNeighbourRanks[i]==neighbourRank){
								exist = 1;
								break;
							}
						}
						if(exist != 1){
							newDomain->listNeighbourRanks[counter] = neighbourRank;
							for(int j=0; j<3; j++) newDomain->neighbourCoords[ndims*counter+j] = (neighbourCoords[j]+newDomain->dims[j])%newDomain->dims[j];
							counter++;
						}
					}
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		for(int j=0; j<ndims; j++)
		{
			tmpDim = 1.f/((float)newDomain->dims[j]);
			newDomain->lowlimit[j] = newDomain->coords[j]*tmpDim;
			newDomain->uplimit[j] = (newDomain->coords[j]+1)*tmpDim;
		}
#endif
	}else{
	  	for(int j=0; j<ndims; j++)
		{
			newDomain->uplimit[j] = 1.;
			newDomain->lowlimit[j] = 0.;
		}
	}
	
	return newDomain;
}

void deallocateDomain(domain_t *thisDomain)
{
	if(thisDomain->listNeighbourRanks != NULL) free(thisDomain->listNeighbourRanks);
	if(thisDomain->neighbourCoords != NULL) free(thisDomain->neighbourCoords);
	free(thisDomain);
}


/*--------------------------------------------------------------------------------------------*/
int getNumNeighbours(int size, int ndims, int *n)
{
	int n2 = n[0];
	int n3 = n[1];
	int numNeighbours = 0;
	if(size == 1)
	{
		numNeighbours = 0;
	}else if(size<27)
	{
		if(size<=12)
		{
			numNeighbours = size-1;//pow(2,n2)+pow(3,n3)-1;
		}else if(size<18){
			size = 11;
		}else numNeighbours = 17;
	}else{
		numNeighbours = pow(3,ndims)-1;
	}
	return numNeighbours;
}

int *getNumCuts(int size)
{
	int n2, n3;
	int *n;
	n2 = 0;
	n3 = 0;
	int tmp = size;
	if(size%2 == 0)
	{
		while((tmp%2 == 0) && (tmp != 0))
		{
			n2++;
			tmp = tmp/2;
		}
		while((tmp%3 == 0) && (tmp != 0))
		{
			n3++;
			tmp = tmp/3;
		}
	}else if(size%3 == 0)
	{
		while((tmp%3 == 0) && (tmp != 0))
		{
			n3++;
			tmp = tmp/3;
		}
		while((tmp%2 == 0) && (tmp != 0))
		{
			n2++;
			tmp = tmp/2;
		}
	}else if(size == 1)
	{
		n2 = 0;
		n3 = 0;
	}
	
	n = malloc(2*sizeof(int));
	n[0] = n2;
	n[1] = n3;
	
	return n;
}

void getDims(int ndims, int *n, int *dims)
{
	int n2 = n[0];
	int n3 = n[1];
	
	while(n2>0 || n3>0)
	{
		for(int i=0; i<ndims; i++)
		{
			if(n3>0)
			{
				dims[i] = dims[i]*3;
				n3--;
			}else if(n2>0)
			{
				dims[i] = dims[i]*2;
				n2--;
			}
		}
	}
}