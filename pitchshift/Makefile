CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: pitchshift

pitchshift: src/pitchshift.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: pitchshift
	./pitchshift

clean:
	rm -f pitchshift
