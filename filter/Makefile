CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: filter

filter: src/filter.c ../fft/src/fft.h ../limiter/src/limiter.h
	$(CC) $(CFLAGS) $< -o $@ -lm

test: filter
	./filter

clean:
	rm -f filter
