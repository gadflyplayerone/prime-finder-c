CC=gcc
CFLAGS=-O3 -march=native -Wall
LDFLAGS= -lgmp

all: prime_finder

prime_finder: main.c
	$(CC) $(CFLAGS) -o prime_finder main.c $(LDFLAGS)

clean:
	rm -f prime_finder