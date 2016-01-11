#ifndef READ_GADGET_FILE_H
#define READ_GADGET_FILE_H
#endif


void read_gadget_file(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles);
void read_particle_pos(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, int particle_type, char buf[]);
void recv_particle_pos(domain_t * thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, int particle_type);
void sort_particles_to_processors(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, float *buf, int Npart_chunk, char description[]);