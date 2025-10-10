CC=gcc
CFLAGS=-O3 -march=native -mtune=native -fno-plt -fopenmp -flto -fomit-frame-pointer -funroll-loops
LDFLAGS= -lgmp -lm -fopenmp -flto

all: prime_finder
prime_finder: main.c
	$(CC) $(CFLAGS) -o prime_finder main.c $(LDFLAGS)

clean:
	rm -f prime_finder
