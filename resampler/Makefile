CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: resampler

resampler: src/resampler.c ../common/sr.h
	$(CC) $(CFLAGS) -I../common $< -o $@ -lm

test: resampler
	./resampler

clean:
	rm -f resampler
