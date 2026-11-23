CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: phasecheck

phasecheck: src/phasecheck.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: phasecheck
	./phasecheck --selftest

clean:
	rm -f phasecheck

.PHONY: all test clean
