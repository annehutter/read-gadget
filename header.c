#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "header.h"

header_t *initHeader()
{
	header_t *newHeader;
	newHeader = malloc(sizeof(header_t));
	if(newHeader == NULL)
	{
		fprintf(stderr, "ERROR: initHeader: Not enough memory to allocate header.\n");
		exit(EXIT_FAILURE);
	}
	
	for(int i=0; i<6; i++)
	{
		newHeader->npart[i] = 0;
		newHeader->mass[i] = 0.;
		newHeader->npartTotal[i] = 0;
	}
	newHeader->time = 0.;
	newHeader->redshift = 0.;
	newHeader->flag_sfr = 0;
	newHeader->flag_feedback = 0;
	newHeader->flag_cooling = 0;
	newHeader->num_files = 0;
	newHeader->BoxSize = 0.;
	newHeader->Omega0 = 0.;
	newHeader->OmegaLambda = 0.;
	newHeader->HubbleParam = 0.;
	
	return newHeader;
}

void deallocateHeader(header_t *thisHeader)
{
	free(thisHeader);
}

void read_header_gadget1(header_t *thisHeader, char *fname, int myRank)
{
  	FILE *fd;
	char buf[200];
	
	int dummy;
	float dummyfloat;

	sprintf(buf, "%s", fname);

	if(!(fd = fopen(buf, "r"))){
			printf("can't open file `%s`\n", buf);
			exit(0);
		}

	printf("header: reading `%s' ...\n", buf);
	fflush(stdout);
	
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	
	fread(thisHeader, sizeof(header_t), 1, fd);
	
	printf("%d\t%d\t%d\t%d\t%d\t%d\n", thisHeader->npart[0], thisHeader->npart[1], thisHeader->npart[2], thisHeader->npart[3], thisHeader->npart[4], thisHeader->npart[5]);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	
	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[1]+thisHeader->npart[4])*3; i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//POS
	}
	printf("rank %d: pos = %e\n", myRank, dummyfloat);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);
	
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);
	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[1]+thisHeader->npart[4])*3; i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	// VEL
	}
	printf("rank %d: vel = %e\n", myRank, dummyfloat);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	
	
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[1]+thisHeader->npart[4]); i++)
	{
		fread(&dummy, sizeof(dummy), 1, fd);
	}
	printf("rank %d: ID = %d\n", myRank, dummy);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	
	fclose(fd);
}

void read_header_gadget2(header_t *thisHeader, char *fname, int myRank)
{
  	FILE *fd;
	char buf[200];
	
	int dummy;
	char dummychar;
	float dummyfloat;
	
	sprintf(buf, "%s", fname);

	if(!(fd = fopen(buf, "r"))){
			printf("can't open file `%s`\n", buf);
			exit(0);
		}

	printf("reading `%s' ...\n", buf);
	fflush(stdout);
	
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf("%c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf("%c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	
	fread(thisHeader, sizeof(header_t), 1, fd);
	printf("%d\t%d\t%d\t%d\t%d\t%d\n", thisHeader->npart[0], thisHeader->npart[1], thisHeader->npart[2], thisHeader->npart[3], thisHeader->npart[4], thisHeader->npart[5]);
	printf("boxsize = %e\n",thisHeader->BoxSize);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: dummy = %d\n",myRank,dummy);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\t%d\n",myRank,sizeof(dummy),dummy, thisHeader->npartTotal[1]);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[1]+thisHeader->npart[4])*3; i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//POS
	}
// 	printf("rank %d: %d dummy = %f\n",myRank,sizeof(dummyfloat),dummyfloat);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[1]+thisHeader->npart[4])*3; i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//VEL
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[1]+thisHeader->npart[4]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[4]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

			fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

		fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

			fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

		fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

		fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

		fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[4]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummychar, 1,1,fd);
	printf("rank %d: %d dummychar = %c",myRank,4,dummychar);
	for(int j=0; j<2; j++){
		fread(&dummychar, 1,1,fd);
		printf(" %c",dummychar);
	}
	fread(&dummychar, 1,1,fd);
	printf(" %c\n",dummychar);
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

		fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);

	for(unsigned int i=0; i<(thisHeader->npart[0]+thisHeader->npart[4]); i++)
	{
		fread(&dummyfloat, sizeof(dummyfloat), 1, fd);	//ID
	}
	fread(&dummy, sizeof(dummy), 1, fd);
	printf("rank %d: %ld dummy = %d\n",myRank,sizeof(dummy),dummy);
	
	fclose(fd);
}