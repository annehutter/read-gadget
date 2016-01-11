#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "part.h"

part_t *initParticles()
{
	part_t *newPart;
	newPart = malloc(sizeof(part_t));
	
	newPart->num = 0;
	newPart->pos = NULL;
	newPart->vel = NULL;
	newPart->ID = NULL;
	newPart->mass = NULL;
	newPart->u = NULL;
	newPart->rho = NULL;
	newPart->ne = NULL;
	newPart->nh = NULL;
	newPart->hsml = NULL;
	newPart->sfr = NULL;
	newPart->age = NULL;
	newPart->z = NULL;
	
	return newPart;
}

void allocateParticles_pos(part_t *theseParticles, int numPart, float *buf)
{
	unsigned int numOld = theseParticles->num;
	theseParticles->num = numOld + numPart;
	
	if(theseParticles->pos == NULL)
	{
		theseParticles->pos = malloc(sizeof(float)*numPart*3);
	}else{
		theseParticles->pos = realloc(theseParticles->pos, sizeof(float)*(theseParticles->num)*3);
	}
	
	memcpy(&theseParticles->pos[3*numOld], buf, sizeof(float)*numPart*3);
// 	for(int p=0; p<3*numPart, p++) theseParticles->pos[3*numOld+p] = buf[p];
}

void deallocateParticles(part_t *theseParticles)
{
	if(theseParticles->pos != NULL) free(theseParticles->pos);
	if(theseParticles->vel != NULL) free(theseParticles->vel);
	if(theseParticles->ID != NULL) free(theseParticles->ID);
	if(theseParticles->mass != NULL) free(theseParticles->mass);
	if(theseParticles->u != NULL) free(theseParticles->u);
	if(theseParticles->rho != NULL) free(theseParticles->rho);
	if(theseParticles->ne != NULL) free(theseParticles->ne);
	if(theseParticles->nh != NULL) free(theseParticles->nh);
	if(theseParticles->hsml != NULL) free(theseParticles->hsml);
	if(theseParticles->sfr != NULL) free(theseParticles->sfr);
	if(theseParticles->age != NULL) free(theseParticles->age);
	if(theseParticles->z != NULL) free(theseParticles->z);
	free(theseParticles);
}