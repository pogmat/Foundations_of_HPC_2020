PGM_DIR = ../../pgm
PGM_OBJ = pgm_bin.o
KERNEL_DIR = ../../kernel
KERNEL_OBJ = kernel.o

CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -march=native -O3
LDFLAGS = -lm

.PHONY: ugly_p
ugly_p: subdirs ugly_p.x

ugly_p.x: ugly_p.o $(PGM_DIR)/$(PGM_OBJ) $(KERNEL_DIR)/$(KERNEL_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

ugly_p.o: ugly_p.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(PGM_DIR)/$(PGM_OBJ):
	$(MAKE) -C $(PGM_DIR)

$(KERNEL_DIR)/$(KERNEL_OBJ):
	$(MAKE) -C $(KERNEL_DIR)

.PHONY: subdirs
subdirs: FORCE
	$(MAKE) -C $(PGM_DIR)
	$(MAKE) -C $(KERNEL_DIR)

.PHONY: FORCE
FORCE:

.PHONY: clean
clean:
	rm -f ugly_p.o ugly_p.x
