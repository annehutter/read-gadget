#ifndef SPH_DENSITY_H
#define SPH_DENSITY_H
#endif

float calc_distance(float x1, float y1, float z1, float x2, float y2, float z2);
float get_SPH_density(float distance, char *kernel);
void compute_SPH_density(domain_t *thisDomain, header_t *thisHeader, part_t *theseParticles);