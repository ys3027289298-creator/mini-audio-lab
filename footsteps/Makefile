CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: footsteps

footsteps: src/footsteps.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: footsteps
	./footsteps --test

# writes footsteps/samples/<material>_<00..09>.wav and samples.csv
samples: footsteps
	mkdir -p samples
	./footsteps --wav samples --csv > samples/samples.csv

clean:
	rm -rf footsteps samples

.PHONY: all test samples clean
