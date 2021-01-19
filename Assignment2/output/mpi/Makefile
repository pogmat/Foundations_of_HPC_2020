 MPICC = mpicc
 CFLAGS = -Wall -Wextra -Wpedantic -O3 -march=native -lm
 
 blur.mpi.x : blur.mpi.c
	$(MPICC) $(CFLAGS) $^ -o $@
	
.PHONY: vignetting
vignetting: blur.mpi.c
	$(MPICC) $(CFLAGS) -DVIGNETTING $^ -o blur.mpi.x
