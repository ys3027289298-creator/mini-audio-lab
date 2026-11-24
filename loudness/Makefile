CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: loudness

loudness: src/loudness.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: loudness
	./loudness

clean:
	rm -f loudness
