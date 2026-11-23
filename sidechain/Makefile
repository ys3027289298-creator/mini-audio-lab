CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: sidechain

sidechain: src/sidechain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sidechain
	./sidechain

clean:
	rm -f sidechain

.PHONY: all test clean
