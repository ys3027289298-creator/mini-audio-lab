CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: stft

stft: src/stft.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: stft
	./stft ../artifacts

clean:
	rm -f stft
