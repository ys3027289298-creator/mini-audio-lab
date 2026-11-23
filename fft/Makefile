CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: fft

fft: src/fft.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: fft
	./fft

clean:
	rm -f fft
