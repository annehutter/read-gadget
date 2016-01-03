#ifndef PART_H
#define PART_H
#endif

typedef struct
{
	int num;
	float *pos;
	float *vel;
	unsigned int *ID;
	float *mass;
	float *u;
	float *rho;
	float *ne;
	float *nh;
	float *hsml;
	float *sfr;
	float *age;
	float *z;
} part_t;

part_t *initParticles();
void allocateParticles_pos(part_t *theseParticles, int numPart, float *buf);
void deallocateParticles(part_t *theseParticles);
