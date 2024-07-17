# Makefile

# Compiler
CC = gcc
MPICC = mpicc

# Compiler flags
CFLAGS = -Wall -O2
MPICFLAGS = -Wall -O2

# Targets
TARGETS = fft seq-dht seqbfly par-dht parbly

# Source files
FFT_SRC = fft.c
SEQ_DHT_SRC = seq-dht.c
SEQBFLY_SRC = seqbfly.c
PAR_DHT_SRC = par-dht.c
PARBFLY_SRC = parbly.c

# Object files
FFT_OBJ = $(FFT_SRC:.c=.o)
SEQ_DHT_OBJ = $(SEQ_DHT_SRC:.c=.o)
SEQBFLY_OBJ = $(SEQBFLY_SRC:.c=.o)
PAR_DHT_OBJ = $(PAR_DHT_SRC:.c=.o)
PARBFLY_OBJ = $(PARBFLY_SRC:.c=.o)

# Rules
all: $(TARGETS)

fft: $(FFT_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

seq-dht: $(SEQ_DHT_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

seqbfly: $(SEQBFLY_OBJ)
	$(CC) $(CFLAGS) -o $@ $^

par-dht: $(PAR_DHT_OBJ)
	$(MPICC) $(MPICFLAGS) -o $@ $^

parbly: $(PARBFLY_OBJ)
	$(MPICC) $(MPICFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

par-dht.o: par-dht.c
	$(MPICC) $(MPICFLAGS) -c $< -o $@

parbly.o: parbly.c
	$(MPICC) $(MPICFLAGS) -c $< -o $@

clean:
	rm -f $(FFT_OBJ) $(SEQ_DHT_OBJ) $(SEQBFLY_OBJ) $(PAR_DHT_OBJ) $(PARBFLY_OBJ) $(TARGETS)

.PHONY: all clean
