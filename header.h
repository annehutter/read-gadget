#ifndef HEADER_H
#define HEADER_H
#endif

typedef struct
{
	unsigned int npart[6];
	double mass[6];
	double time;
	double redshift;
	int flag_sfr;
	int flag_feedback;
	unsigned int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
} header_t;

header_t *initHeader();
void read_header_gadget1(header_t *thisHeader, char *fname, int myRank);
void read_header_gadget2(header_t *thisHeader, char *fname, int myRank);
void deallocateHeader(header_t *thisHeader);