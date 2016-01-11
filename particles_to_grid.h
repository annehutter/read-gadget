#ifndef PARTICLES_TO_GRID_H
#define PARTICLES_TO_GRID_H
#endif

double meanRho_part(part_t *theseParticles);
double meanRho_part_global(part_t *theseParticles);
void comp_arrays_for_clump(domain_t *thisDomain, header_t *thisHeader, part_t *theseParticles, grid_t *thisGrid);
void produce_clumping_factor_fields(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, int gridsize);