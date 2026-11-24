CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -I../stream

all: resampler

resampler: src/resampler.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: resampler
	./resampler

clean:
	rm -f resampler
