PGM_DIR = ../pgm
PGM_OBJ = pgm_bin.o
KERNEL_DIR = ../kernel
KERNEL_OBJ = kernel.o

CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic
LDFLAGS =

.PHONY: ugly
ugly: subdirs ugly.x

ugly.x: ugly.o $(PGM_DIR)/$(PGM_OBJ) $(KERNEL_DIR)/$(KERNEL_OBJ)
	$(CC) $(LDFLAGS) $^ -o $@

ugly.o: ugly.c
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
	rm -f ugly.o ugly.x
