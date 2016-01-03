#ifndef INPUT_H
#define INPUT_H
#endif

typedef struct
{
	int snapshot_format;
	int files;
	char *fname;
	
	int particle_type;
	
	int size_in_MB;
} input_t;


input_t *initInput(int files, char *fname, int snapshot_format, int particle_type, int size_in_MB);
void deallocateInput(input_t *thisInput);
