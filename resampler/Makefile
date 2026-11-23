CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: resampler aliasing_bound

resampler: src/resampler.c
	$(CC) $(CFLAGS) $< -o $@ -lm

aliasing_bound: src/aliasing_bound.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: resampler aliasing_bound
	./resampler
	./aliasing_bound

clean:
	rm -f resampler aliasing_bound
