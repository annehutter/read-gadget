#ifndef GRID_GLOBAL_H
#define GRID_GLOBAL_H
#endif


grid_t *initGrid_withDomain(domain_t *thisDomain, float boxsize, int gridsize);
float *allocateGrid_withDomain_float(domain_t *thisDomain, int gridsize);
int *allocateGrid_withDomain_int(domain_t *thisDomain, int gridsize);

#ifdef __MPI
void save_rho_to_file(domain_t *thisDomain, grid_t *thisGrid, input_t *thisInput);
void save_inv_rho_to_file(grid_t *thisGrid, input_t *thisInput);
void save_npart_cell_to_file(grid_t *thisGrid, input_t *thisInput);
#endif