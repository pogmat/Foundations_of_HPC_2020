SOURCE = kernel.c
OBJECT = $(SOURCE:.c=.o)
HEADER = $(SOURCE:.c=.h)

TESTSRC = test.c
TESTOBJ = $(TESTSRC:.c=.o)
TESTTRG = $(TESTSRC:.c=.x)

CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -march=native -O0 -g
LDFLAGS =

$(OBJECT): $(SOURCE) $(HEADER)
	$(CC) -c $(CFLAGS) -o $@ $<

$(TESTOBJ): $(TESTSRC) $(HEADER)
	$(CC) -c $(CFLAGS) -o $@ $<

$(TESTTRG): $(OBJECT) $(TESTOBJ)
	$(CC) $(LDFLAGS) $^ -o $@

.PHONY: test
test: $(TESTTRG)

.PHONY: clean
clean:
	rm -f $(OBJECT) $(TESTOBJ) $(TESTTRG)
