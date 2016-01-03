#ifndef DOMAIN_DECOMPOSITION_H
#define DOMAIN_DECOMPOSITION_H
#endif


typedef struct
{
	int ndims;
	int dims[3];
	int periods[3];
	int coords[3];
	
	int originRank;
	int size;

	int numNeighbours;
	int *listNeighbourRanks;
	int *neighbourCoords;
	
	float uplimit[3];
	float lowlimit[3];
} domain_t;

domain_t *initDomain(int thisRank, int size);
void deallocateDomain(domain_t *thisDomain);

int getNumNeighbours(int size, int ndims, int *n);
int *getNumCuts(int size);
void getDims(int ndims, int *n, int *dims);
