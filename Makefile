CC=gcc
#CFLAGS=-O3 -march=native -Wall
CFLAGS=-O3 -march=native -mtune=native -fno-plt #x86
#CFLAGS=-O3 -mcpu=native #graviton arm64
# Add -lm to LDFLAGS
LDFLAGS= -lgmp -lm

all: prime_finder

prime_finder: main.c
	$(CC) $(CFLAGS) -o prime_finder main.c $(LDFLAGS)

clean:
	rm -f prime_finder