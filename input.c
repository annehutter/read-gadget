#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "input.h"

input_t *initInput(int files, char *fname, int snapshot_format, int particle_type, int size_in_MB)
{
	input_t *newInput;
	newInput = malloc(sizeof(input_t));
		
	newInput->files = files;
	newInput->fname = fname;
	newInput->snapshot_format = snapshot_format;
	newInput->particle_type = particle_type;
	
	newInput->size_in_MB = size_in_MB;
	
	return newInput;
}

void deallocateInput(input_t *thisInput)
{
	free(thisInput);
}