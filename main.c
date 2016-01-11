#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "input.h"
#include "header.h"
#include "domain_decomposition.h"
#include "part.h"
#include "grid.h"

#include "read_gadget_file.h"
#include "sph_density.h"
#include "particles_to_grid.h"

int main(int argc, char **argv)
{
	int size=1;
	int myRank=0;
	
	input_t *input;
	
	header_t *header;
	header = initHeader();
	
	part_t *particles;
	particles = initParticles();

	int snapshot_format;

#ifdef __MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

	domain_t *domain;
	domain = initDomain(myRank, size);
		
	/* get snapshotfilenames */
	if(argc!=8 && myRank==0){
		printf("usage: ./read_snapshot_parallel <path of snapshot> <basename of snapshot> <snapshot number> <number of files per snapshot> <size in MB of chunks> <gridsize of clump array>\n");
		exit(0);
	}
	char path[200], input_fname[200], basename[200];
	int snapshot_number, files, size_in_MB;
	int gridsize;

	sprintf(path, argv[1]);
	sprintf(basename, argv[2]);
	snapshot_number = atoi(argv[3]);
	files = atoi(argv[4]);
	snapshot_format = atoi(argv[5]);
	size_in_MB = atoi(argv[6]);
	gridsize = atoi(argv[7]);
	
	if(gridsize%size != 0)
	{
		if(myRank == 0) printf("number of processors needs to be a multiple of gridsize (i.e. 2 and 3) are allowed. Please adjust your number of processors or gridsize!\n");
		exit(0);
	}

	sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
	
	input = initInput(files, input_fname, snapshot_format, 1, size_in_MB);
	
	/* read particles */
	
	read_gadget_file(domain, header, input, particles);
	
	/* compute SPH density for each particle */
	compute_SPH_density(domain, header, particles);
	
	produce_clumping_factor_fields(domain, header, input, particles, gridsize);

	printf("deallocating ...\n");
	deallocateDomain(domain);
	printf("- domain\n");
	deallocateHeader(header);
	printf("- header\n");
	deallocateParticles(particles);
	printf("- particles\n");
	deallocateInput(input);
	printf("- input\n");

#ifdef __MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	
	return 0;
}