CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

BINS = fft spectral

all: $(BINS)

fft: src/fft.c src/fft_core.c src/fft_internal.h
	$(CC) $(CFLAGS) src/fft.c src/fft_core.c -o $@ -lm

spectral: src/spectral.c src/fft_core.c src/fft_internal.h
	$(CC) $(CFLAGS) src/spectral.c src/fft_core.c -o $@ -lm

test: $(BINS)
	./fft
	./spectral

clean:
	rm -f $(BINS)
