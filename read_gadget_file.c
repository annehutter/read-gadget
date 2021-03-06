#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "input.h"
#include "header.h"
#include "domain_decomposition.h"
#include "part.h"
#include "read_gadget_file.h"

void read_gadget_file(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles)
{
	char buf[200];
	
	int myRank = thisDomain->originRank;
	int size = thisDomain->size;
	
	int files = thisInput->files;
	char *fname = thisInput->fname;
	int snapshot_format = thisInput->snapshot_format;
	int particle_type = thisInput->particle_type;
	
	int max;
	
	if(size<files)
	{
		if(files%size == 0) max = files;
		else max = files + (size - files%size);
	}
	else max = size;
	
	for(int i=myRank%size; i<max; i=i+size)
	{
		if(i < files)
		{
			if(files > 1)
			{
				sprintf(buf, "%s.%d", fname, i);
			}else{
				sprintf(buf, "%s", fname);
			}

			
	/*-----------------------------------------------------------------------------------------*/
	/* reading header of file                                                                  */
	/*-----------------------------------------------------------------------------------------*/

			if(snapshot_format == 1)
			{
				read_header_gadget1(thisHeader, buf, myRank);
			}else if(snapshot_format == 2)
			{
				read_header_gadget2(thisHeader, buf, myRank);
			}else{
				fprintf(stderr,"Given snapshot format not supported.\n");
			}
		}

#ifdef __MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

/*-----------------------------------------------------------------------------------------*/
/* reading positions of particles (in data chunks)                                         */
/*-----------------------------------------------------------------------------------------*/
		if(i<files)
		{
			read_particle_pos(thisDomain, thisHeader, thisInput, theseParticles, particle_type, buf);
		}else{
			recv_particle_pos(thisDomain, thisHeader, thisInput, theseParticles, particle_type);
		}
	}
}


void read_particle_pos(domain_t * thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, int particle_type, char buf[])
{
	int file_offset;
	
	unsigned int buf_Npos;
	float *buf_pos;
	char buf_char[4];
	unsigned int Npart_chunk;
	unsigned int Nchunks=0, maxNchunks=0;
	int allParticleRead = 0;
	int particleRead = 0;
	
	FILE *fd;
	
	unsigned int numTotPart;
	
	int snapshot_format = thisInput->snapshot_format;
	
	if(!(fd = fopen(buf, "r"))){
		printf("can't open file `%s`\n", buf);
		exit(0);
	}

	printf("reading `%s' ...\n", buf);
	fflush(stdout);

	numTotPart = 0;
	for(int i=0; i<6; i++) 
	{
		numTotPart += thisHeader->npart[i];
		printf("type %d: %d\n", i, thisHeader->npart[i]);
	}
	
	/* skip header */
	if(snapshot_format == 1)
	{
		file_offset = 4+256+4;
		fseek(fd, file_offset, SEEK_SET);

	}else if(snapshot_format == 2)
	{
		file_offset = 4*5+256+4 + 4;
		fseek(fd, file_offset, SEEK_SET);
		
		fread(&buf_char, sizeof(char), 4, fd);
		assert(strcmp(buf_char, "POS ") == 0);
		file_offset = 2*4;
		fseek(fd, file_offset, SEEK_CUR);
	}
	
	/* read number of bytes that hold positions */
	fread(&buf_Npos, sizeof(int), 1, fd);
	assert(buf_Npos == numTotPart*sizeof(float)*3);
	
	/* seek position of particle type in file */
	file_offset = 0;
	for(int i=0; i<particle_type; i++) file_offset += thisHeader->npart[i]*sizeof(float);
	fseek(fd, file_offset, SEEK_CUR);
	
	Npart_chunk = thisInput->size_in_MB*1024*1024/12;
	printf("Npart_chunk = %d corrsponds to %ld bytes\n", Npart_chunk, Npart_chunk*3*sizeof(float));
	buf_pos = malloc(Npart_chunk*3*sizeof(float));
	
	Nchunks = thisHeader->npart[particle_type]/Npart_chunk+1;
#ifdef __MPI
	MPI_Allreduce(&Nchunks, &maxNchunks, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
	maxNchunks = Nchunks;
#endif
	printf("rank %d: Nchunks = %d \tmaxNchunks = %d\n", thisDomain->originRank, Nchunks, maxNchunks);
	allParticleRead = 0;

	for(unsigned int i=0; i<maxNchunks; i++)
	{
		if(allParticleRead == 1) Npart_chunk = 0;
		if((i+1)*Npart_chunk >= thisHeader->npart[particle_type])
		{
			Npart_chunk = thisHeader->npart[particle_type] - i*Npart_chunk;
			buf_pos = realloc(buf_pos, Npart_chunk*3*sizeof(float));
			allParticleRead = 1;
		}
		particleRead += Npart_chunk;

		/* read data in chunks */
		printf("rank %d: reading %d. chunk ... < %d \n", thisDomain->originRank,particleRead, thisHeader->npart[particle_type]);
		fread(buf_pos, sizeof(float), Npart_chunk*3, fd);

#ifdef __MPI
		/* sort particles to processors */
		printf("sorting particles to processores ...\n");
		sort_particles_to_processors(thisDomain, thisHeader, thisInput, theseParticles, buf_pos, Npart_chunk, "POS ");
#else
		allocateParticles_pos(theseParticles, Npart_chunk, buf_pos);
#endif
	}
	free(buf_pos);
	buf_pos = NULL;

	fclose(fd);
	printf("closed file\n");
}

void recv_particle_pos(domain_t * thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, int particle_type)
{
	float *buf_pos;
	unsigned int Npart_chunk;
	unsigned int Nchunks=0, maxNchunks=0;
	int allParticleRead = 0;
	
	Npart_chunk = 0;
	buf_pos = malloc(Npart_chunk*3*sizeof(float));
	
	Nchunks = 0;
#ifdef __MPI
	MPI_Allreduce(&Nchunks, &maxNchunks, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
	maxNchunks = Nchunks;
#endif
	printf("rank %d: Nchunks = %d \tmaxNchunks = %d\n", thisDomain->originRank, Nchunks, maxNchunks);
	allParticleRead = 0;

	for(unsigned int i=0; i<maxNchunks; i++)
	{
#ifdef __MPI
		/* sort particles to processors */
		printf("sorting particles to processores ...\n");
		sort_particles_to_processors(thisDomain, thisHeader, thisInput, theseParticles, buf_pos, Npart_chunk, "POS ");
#else
		allocateParticles_pos(theseParticles, Npart_chunk, buf_pos);
#endif
	}
	free(buf_pos);
	buf_pos = NULL;
}

/* routine to sort read particle position data to processors depending on their position */
void sort_particles_to_processors(domain_t *thisDomain, header_t *thisHeader, input_t *thisInput, part_t *theseParticles, float *buf, int Npart_chunk, char description[])
{
#ifdef __MPI
	int *dims = thisDomain->dims;
	float boxsize = thisHeader->BoxSize;
	float inv_boxsize = 1.f/boxsize;
	
	int *buf_proc;
	int *NpartDomain, *recNpartDomain;
	
	float xpos, ypos, zpos;
	int x, y, z;
	
	float *buf_part_domain, *rec_buf_part_domain;
	
	int counter;
	
	MPI_Status status;
	
	NpartDomain = malloc(thisDomain->size*sizeof(int));
	if(NpartDomain == NULL)
	{
		fprintf(stderr, "ERROR: allocating NpartDomain: Not enough memory to allocate array.\n");
		exit(EXIT_FAILURE);
	}
	
	recNpartDomain = malloc(thisDomain->size*sizeof(int));
	if(recNpartDomain == NULL)
	{
		fprintf(stderr, "ERROR: allocating recNpartDomain: Not enough memory to allocate array.\n");
		exit(EXIT_FAILURE);
	}
	
	buf_proc = malloc(Npart_chunk*sizeof(float));
	if(buf_proc == NULL)
	{
		fprintf(stderr, "ERROR: allocating buf_proc: Not enough memory to allocate array.\n");
		exit(EXIT_FAILURE);
	}
	
	
	for(int i=0; i<thisDomain->size; i++)
	{
		NpartDomain[i] = 0;
	}
	
	for(int p=0; p<Npart_chunk; p++)
	{
		xpos = buf[p*3];
		ypos = buf[p*3+1];
		zpos = buf[p*3+2];

		x = floor(xpos*inv_boxsize*dims[0]);
		y = floor(ypos*inv_boxsize*dims[1]);
		z = floor(zpos*inv_boxsize*dims[2]);
		
		if(x>=thisDomain->size || y>=thisDomain->size || z>=thisDomain->size) printf("P %d: x=%d\t y=%d\t z=%d\t %f %f %f\n",p,x,y,z,xpos,ypos,zpos);
		buf_proc[p] = x*dims[1]*dims[2]+y*dims[2]+z;
		assert(buf_proc[p] == x*dims[1]*dims[2]+y*dims[2]+z);
		assert(buf_proc[p]<thisDomain->size);

		NpartDomain[buf_proc[p]]++;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* sending particles to each domain from originRank */
	for(int domain=0; domain<thisDomain->size; domain++)
	{
		buf_part_domain = malloc(NpartDomain[domain]*3*sizeof(float));
		counter = 0;
		for(int p=0; p<Npart_chunk; p++)
		{
			if(buf_proc[p] == domain)
			{
				buf_part_domain[counter*3] = buf[p*3];
				buf_part_domain[counter*3+1] = buf[p*3+1];
				buf_part_domain[counter*3+2] = buf[p*3+2];
				counter++;
			}
		}
			
		if(domain != thisDomain->originRank)
		{
			MPI_Ssend(&NpartDomain[domain], 1, MPI_INT, domain, 100, MPI_COMM_WORLD);
			MPI_Ssend(buf_part_domain, 3*NpartDomain[domain], MPI_FLOAT, domain, 101, MPI_COMM_WORLD);
		}else{
			if(strcmp(description, "POS ") == 0) allocateParticles_pos(theseParticles, NpartDomain[domain], buf_part_domain);
		}
			
		free(buf_part_domain);
		/* receving particles from domains to OriginRank */
			
		if(domain == thisDomain->originRank){
			for(int recDomain=0; recDomain<thisDomain->size; recDomain++)
			{
				if(recDomain != thisDomain->originRank)
				{
					MPI_Recv(&recNpartDomain[recDomain],1,MPI_INT,recDomain,100,MPI_COMM_WORLD,&status);
					rec_buf_part_domain = malloc(recNpartDomain[recDomain]*3*sizeof(float));
					MPI_Recv(rec_buf_part_domain,3*recNpartDomain[recDomain],MPI_FLOAT,recDomain,101,MPI_COMM_WORLD,&status);
					
					if(strcmp(description, "POS ") == 0) allocateParticles_pos(theseParticles, recNpartDomain[recDomain], rec_buf_part_domain);
					
					free(rec_buf_part_domain);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	free(NpartDomain);
	free(recNpartDomain);
	free(buf_proc);
#endif
}