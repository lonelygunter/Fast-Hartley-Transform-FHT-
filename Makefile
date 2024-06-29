# Define the compiler and flags
CC = gcc
CFLAGS = -Wall -O2

# Define the source and object files
SRC = fft.c seq-dht.c
OBJS = $(SRC:.c=.o)

# Default target
all: fft seq-dht

# Rule to build the target executable
$@: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS)

# Rule to build the object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to clean the build files
clean:
	rm -f $@ $(OBJS)

# Phony targets
.PHONY: all clean