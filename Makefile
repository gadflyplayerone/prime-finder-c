UNAME_S := $(shell uname -s)

CC ?= gcc
CC_NOOMP ?= cc

SRC = main.c
LIBS = -lgmp -lm
OMPFLAGS ?= -fopenmp

CFLAGS_WARN := -Wall -Wextra
CFLAGS_OPT := -O3 -fomit-frame-pointer -funroll-loops
CFLAGS_NATIVE := $(CFLAGS_WARN) $(CFLAGS_OPT) -flto

ifeq ($(UNAME_S),Darwin)
GMP_PREFIX ?= /opt/homebrew
CFLAGS_NOOMP := $(CFLAGS_WARN) $(CFLAGS_OPT)
LDFLAGS_NOOMP :=
GMP_INCLUDE_DIR := $(GMP_PREFIX)/include
GMP_LIB_DIR := $(GMP_PREFIX)/lib
ifneq ($(wildcard $(GMP_INCLUDE_DIR)),)
CFLAGS_NOOMP += -I$(GMP_INCLUDE_DIR)
endif
ifneq ($(wildcard $(GMP_LIB_DIR)),)
LDFLAGS_NOOMP += -L$(GMP_LIB_DIR)
endif
else
CFLAGS_NATIVE += -march=native -mtune=native -fno-plt
CFLAGS_NOOMP := $(CFLAGS_WARN) $(CFLAGS_OPT) -march=native -mtune=native
endif

CFLAGS_OMP := $(CFLAGS_NATIVE)

all: prime_finder

prime_finder: $(SRC)
	$(CC) $(CFLAGS_NATIVE) -o $@ $(SRC) $(LIBS)

prime_finder_noomp: $(SRC)
	$(CC_NOOMP) $(CFLAGS_NOOMP) -o $@ $(SRC) $(LIBS) $(LDFLAGS_NOOMP)

prime_finder_omp8: $(SRC)
	$(CC) $(CFLAGS_OMP) $(OMPFLAGS) -DDEFAULT_OMP_THREADS=8 -o $@ $(SRC) $(LIBS) $(OMPFLAGS)

prime_finder_omp16: $(SRC)
	$(CC) $(CFLAGS_OMP) $(OMPFLAGS) -DDEFAULT_OMP_THREADS=16 -o $@ $(SRC) $(LIBS) $(OMPFLAGS)

prime_finder_omp32: $(SRC)
	$(CC) $(CFLAGS_OMP) $(OMPFLAGS) -DDEFAULT_OMP_THREADS=32 -o $@ $(SRC) $(LIBS) $(OMPFLAGS)

prime_finder_omp64: $(SRC)
	$(CC) $(CFLAGS_OMP) $(OMPFLAGS) -DDEFAULT_OMP_THREADS=64 -o $@ $(SRC) $(LIBS) $(OMPFLAGS)

prime_finder_omp192: $(SRC)
	$(CC) $(CFLAGS_OMP) $(OMPFLAGS) -DDEFAULT_OMP_THREADS=192 -DPAR_CHUNK_SIZE=256 -o $@ $(SRC) $(LIBS) $(OMPFLAGS)

clean:
	rm -f prime_finder prime_finder_noomp prime_finder_omp8 prime_finder_omp16 prime_finder_omp32 prime_finder_omp64 prime_finder_omp192 heatmap_exclusions.ppm heatmap_gcd.ppm heatmap.csv heatmap.ppm results.txt batch-results.txt
