#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include <kdtree.h>

#include "input.h"
#include "header.h"
#include "domain_decomposition.h"
#include "part.h"

#include "sph_density.h"

#define INV_PI 0.3183098861837907

float calc_distance(float x1, float y1, float z1, float x2, float y2, float z2)
{
	float tmp_x, tmp_y, tmp_z;
	
	tmp_x = x1-x2;
	tmp_y = y1-y2;
	tmp_z = z1-z2;
	
	return sqrt(tmp_x*tmp_x+tmp_y*tmp_y+tmp_z*tmp_z);
}

float get_SPH_density(float r, float h, char *kernel)
{
	/* distance is given in units of the smoothing length h */
	/* result is density in units of h^-3 */
	
	float distance = r/h;
	float tmp, factor;
	
	if(strcmp(kernel, "CIC") == 0)
	{
		if(distance>=1) tmp = 0.;
		else tmp = 0.375*INV_PI;
	}
	else if(strcmp(kernel, "TSC") ==0)
	{
		if(distance>=1) tmp = 0.;
		else tmp = 3.*INV_PI*(1.-distance);
	}
	else if(strcmp(kernel, "SPH") == 0)
	{
		if(distance > 1.) 
		{
			tmp = 0.;
// 			printf("tmp2 = %e\n",tmp);
		}else{
			factor = INV_PI*8./(h*h*h);
			if(distance >= 0.5) tmp = factor*(1.-6.*distance*distance*(1.-distance));
			else
			{
				tmp = factor*2.*(1.-distance)*(1.-distance)*(1.-distance);
			}
// 			printf("tmp = %e\n",tmp);
		}
	}else tmp = 0.;
	
	return tmp;
}

static float dist_sq( float *a1, float *a2, int dims ) {
  float dist_sq = 0, diff;
  while( --dims >= 0 ) {
    diff = (a1[dims] - a2[dims]);
    dist_sq += diff*diff;
  }
  return dist_sq;
}

void compute_SPH_density(domain_t *thisDomain, header_t *thisHeader, part_t *theseParticles)
{
	float distance, dist_sq;
	float h;
	
	theseParticles->rho = malloc(sizeof(float)*theseParticles->num);
	
	void *kd = kd_create(3);
	
	for(int p=0; p<theseParticles->num; p++)
	{
		kd_insertf(kd, &theseParticles->pos[3*p], 0);
	}
	
	printf("rank %d: done buliding tree\n", thisDomain->originRank);
	
	struct kdres *res;
	char *pch;
	float pos[3], dist;
	float pt[3] = { 0, 0, 1 };
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for(int p=0; p<theseParticles->num; p++)
	{
		res = kd_nearest_nf(kd, &theseParticles->pos[3*p], 40);
		if(p%100000==0 || p<100) printf("rank %d: p= %d\t res = %d\n",thisDomain->originRank,p,kd_res_size(res));
		theseParticles->rho[p] = 0.;
// 		printf("size = %i\n", kd_res_size(res));
		h = 10.;//sqrt(kd_res_item_dist_sq(res, 39));
		for(int i=0; i<40; i++)
		{
			dist_sq = kd_res_item_dist_sq(res, i);
			distance = sqrt(dist_sq);
// 			printf("%i: %e\t %e\t%e\n", i, distance, h,theseParticles->rho[p]);
			theseParticles->rho[p] += get_SPH_density(distance, h, "SPH");
		}
		kd_res_free(res);
		if(theseParticles->rho[p]<=0.) printf("rank %d: p= %d\t rho = %e\n",thisDomain->originRank,p,theseParticles->rho[p]);
		if(p%100000==0 || p<100) printf("rank %d: p= %d\t rho = %e\t %e\t%e\t%e\n",thisDomain->originRank,p,theseParticles->rho[p], theseParticles->pos[3*p], theseParticles->pos[3*p+1], theseParticles->pos[3*p+2]);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	kd_free(kd);
}