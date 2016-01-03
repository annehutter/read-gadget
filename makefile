CC = mpicc
CFLAGS = -c -std=c99 -D __MPI -march=native -O3 -Wall -Wextra -Wshadow -ftree-vectorize 
LDFLAGS = -lm -lmpich
SOURCES = 	./main.c \
		./header.c \
		./read_gadget_file.c \
		./domain_decomposition.c \
		./input.c \
		./part.c
OBJECTS = $(SOURCES:.c=.o)
DOBJECTS = $(SOURCES:.c=.d)
EXECUTABLE = compSPHdens

.PHONY: all clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
