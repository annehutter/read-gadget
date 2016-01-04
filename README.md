# read-gadget
reading gadget files and distributing particles on multiple processors &amp; computing SPH densities and clumping factor

Requirements:
- Serial run: -
- Parallel run: mpi libraries

Compilation: make

Usage: 
- Serial: ./compSPHdens "PATH to directory where files are located" "basename of Gadget file(s)" "snapshot number" "number of files in one snapshot" "Gadget output type (1 or 2)" "size of reading chunks in MB"
- Parallel: mpiexec -np "numProcessors" ./compSPHdens "PATH to directory where files are located" "basename of Gadget file(s)" "snapshot number" "number of files in one snapshot" "Gadget output type (1 or 2)" "size of reading chunks in MB"

