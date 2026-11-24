CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
BLOCKING ?= 1

.PHONY: all bench test clean

all: fft

fft: src/fft.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# block-streaming analyzer benchmark; BLOCKING=0 disables blocking and
# routes every run through the full-segment FFT path
bench: bench/bench_blockfft.c src/blockfft.c
	$(CC) $(CFLAGS) -D_POSIX_C_SOURCE=199309L -DBLOCKING=$(BLOCKING) -Isrc $< -o bench/bench_blockfft -lm

test: fft
	./fft

clean:
	rm -f fft bench/bench_blockfft
