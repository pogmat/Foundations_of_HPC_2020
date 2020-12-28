TARGET = omp_blur

HEADS = grid kernel pgm_bin
SRC = omp_blur grid kernel pgm_bin
LIBS = -lm

TDIR = bin
IDIR = include
BDIR = build
SDIR = src

HEADEXT = h
SRCEXT = c
OBJEXT = o
TRGEXT = x

COMPILEFLAGS =
RUNOPTIONS =

DEPS = $(patsubst %, $(IDIR)/%.$(HEADEXT), $(HEADS))
OBJS = $(patsubst %, $(BDIR)/%.$(OBJEXT), $(SRC))
SCRS = $(patsubst %, $(SDIR)/%.$(SRCEXT), $(SRC))
RULES = $(patsubst %, $(BDIR)/%.d, $(SRC))

ifeq ($(strip $(TRGEXT)), )
	TRGNAME=$(TARGET)
else
	TRGNAME=$(TARGET).$(TRGEXT)
endif

CC = gcc
CFLAGS = -fopenmp -flto -fuse-linker-plugin -fivopts -Wall -Wextra -Wpedantic
CPPFLAGS = $(COMPILEFLAGS)
CPPFLAGS += -I $(IDIR)
LDFLAGS = $(LIBS) -fopenmp
ifeq (${DEBUG}, y)  
	CPPFLAGS += -DDEBUG_MODE
	CFLAGS += -g -O0 -DDEBUG_MODE
else
	CFLAGS += -O3 -march=native
endif

CPPFLAGS += -MMD -MP

.SUFFIXES :

.PHONY : all
all: $(TDIR)/$(TRGNAME)

$(BDIR)/%.$(OBJEXT) : $(SDIR)/%.$(SRCEXT)
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

$(TDIR)/$(TRGNAME): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^

.PHONY : clean
clean :
	rm -f $(TDIR)/$(TRGNAME) $(OBJS) $(RULES)

.PHONY : run
run : $(TDIR)/$(TRGNAME)
	@$(TDIR)/$(TRGNAME) $(RUNOPTIONS)

.PHONY : debug
debug :
	@$(MAKE) DEBUG=y --no-print-directory -B

.PHONY : valgrind
valgrind : $(TDIR)/$(TRGNAME)
	valgrind --leak-check=full --track-origins=yes $(TDIR)/$(TRGNAME) $(RUNOPTIONS)

$(shell mkdir -p $(BDIR) $(TDIR))

-include $(RULES)
