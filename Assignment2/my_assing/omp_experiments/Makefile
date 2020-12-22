CC = gcc
CFLAGS = -fopenmp -Wall -Wextra -Wpedantic -march=native -O3
LDFLAGS = -lm

naive.x: naive.c
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@
