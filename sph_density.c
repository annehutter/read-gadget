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

float get_SPH_density(float distance, char *kernel)
{
	/* distance is given in units of the smoothing length h */
	/* result is density in units of h^-3 */
	
	float tmp, tmp1, tmp2;
	
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
		if(distance >= 2) tmp = 0.;
		else
		{
			tmp1 = 0.25*(2.-distance)*(2.-distance)*(2.-distance);
			if(distance >= 1) tmp = tmp1*INV_PI;
			else
			{
				tmp2 = (1.-distance)*(1.-distance)*(1.-distance);
				tmp = (tmp1 - tmp2)*INV_PI;
			}
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
// 	float xp, yp, zp;
// 	float xp2, yp2, zp2;
// 	
// 	float distance = 0.;
// 	float density;
// 	float *rho;
	
// 	int *gridcells;
// 	int nbins = 128;

	
// 	float boxsize = thisHeader->BoxSize;
// 	float inv_boxsize = 1.f/boxsize;
	
	theseParticles->rho = malloc(sizeof(float)*theseParticles->num);
// 	gridcells = malloc(sizeof(int)*theseParticles->num);	
	
	void *kd = kd_create(3);
// 	
// 	for(int p=0; p<theseParticles->num; p++)
// 	{
// 		kd_insertf(kd, &theseParticles->pos[p], 0);
// 	}
	
	printf("rank %d: done buliding tree\n", thisDomain->originRank);
	
	struct kdres *res;
	char *pch;
	float pos[3], dist;
	float pt[3] = { 0, 0, 1 };
	for(int p=0; p<theseParticles->num; p++)
	{
// 		res = kd_nearest_rangef(kd, &theseParticles->pos[p], 0.1);//, 40);
// 		if(p%100000==0) printf("rank %d: p= %d\t res = %d\n",thisDomain->originRank,p,kd_res_size(res));
		theseParticles->rho[p] = 1.;//kd_res_size(res);
// 		kd_res_free(res);
	}
// 	printf("found %d results:\n", kd_res_size(res));
// 	while( !kd_res_end( res ) ) 
// 	{
// 		/* get the data and position of the current result item */
// 		pch = (char*)kd_res_itemf( res, pos );
// 
// 		/* compute the distance of the current result from the pt */
// 		dist = sqrt( dist_sq( &theseParticles->pos[0], pos, 3 ) );
// 
// 		/* print out the retrieved data */
// 		printf( "node at (%.3f, %.3f, %.3f) is %.3f away\n", pos[0], pos[1], pos[2], dist);
// 
// 		/* go to the next entry */
// 		kd_res_next( res );
// 	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	kd_free(kd);
// 	for(int p=0; p<theseParticles->num; p++)
// 	{
// 		if(p%1000 == 0) printf("rank %d: p = %d\t %f\n", thisDomain->originRank, p/1000, density);
// 		rho[p] = 0.;
// 		
// 		xp = theseParticles->pos[3*p];
// 		yp = theseParticles->pos[3*p+1];
// 		zp = theseParticles->pos[3*p+2];
// 		
// 		gridcells[p] = xp*inv_boxsize*nbins*nbins+yp*inv_boxsize*nbins+zp*inv_boxsize;
// 		density = 0.;
// 		for(int p2=0; p2<theseParticles->num; p2++)
// 		{
// 			if(gridcells[p] == gridcells[p2])
// 			{
// 				xp2 = theseParticles->pos[3*p2];
// 				yp2 = theseParticles->pos[3*p2+1];
// 				zp2 = theseParticles->pos[3*p2+2];
// 				
// 				distance = xp2;
// // 				distance = calc_distance(xp, yp, zp, xp2, yp2, zp2);
// 				density += get_SPH_density(distance, "SPH");
// 			}
// 		}
// // 		rho[p] = density;
// 	}
	
	
		
// 	printf("start sorting\n");
// 	for(int p=1; p<theseParticles->num; p++)
// 	{
// 		if(p%100000==0) printf("%d\n",p);
// 		valueToSort = gridcells[sortedArray[p]];
// 		valueToSort_index = sortedArray[p];
// 		int p2 = p;
// 		while(p2>0 && gridcells[sortedArray[p2-1]]>valueToSort)
// 		{
// 			sortedArray[p2] = sortedArray[p2-1];
// 			p2 = p2-1;
// 		}
// 		sortedArray[p2] = valueToSort_index;
// 	}
// 	printf("end sorting\n");
	
// 	for(int p=0; p<theseParticles->num; p++) if(p<100 && thisDomain->originRank==0) printf("%d\n",sortedArray[p]);

// 		range = 1;
// 		for(int x=-range; x<=range; x++)
// 		{
// 			for(int y=-range; y<=range; y++)
// 			{
// 				for(int z=-range; z<=range; z++)
// 				{
// 					neighbour_cell = x*nbins*nbins+y*nbins+z;
// 					for(int p2=0; p2<theseParticles->num; p2++)
// 					{
// 						if(p2 == neighbour_cell)
// 						{
// 							xp2 = theseParticles->pos[3*p2];
// 							yp2 = theseParticles->pos[3*p2+1];
// 							zp2 = theseParticles->pos[3*p2+2];
// 							
// 							distance = calc_distance(xp, yp, zp, xp2, yp2, zp2);
// 							density = get_SPH_density(distance, "SPH");
// 							
// 							theseParticles->rho[p] += density;
// 						}
// 					}
// 				}
// 			}
// 		}
// 		for(int p2=0; p2<theseParticles->num; p2++)
// 		{
// 			if(p != p2)
// 			{
// 				xp2 = theseParticles->pos[3*p2];
// 				yp2 = theseParticles->pos[3*p2+1];
// 				zp2 = theseParticles->pos[3*p2+2];
// 				
// 				distance = calc_distance(xp, yp, zp, xp2, yp2, zp2);
// 				if(distance<1.)
// 				{
// 					density = get_SPH_density(distance, "SPH");
// 				
// 					theseParticles->rho[p] += density;
// 				}
// 			}
// 		}

// 	free(rho);
// 	free(gridcells);
}