CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: fxchain

fxchain: src/fxchain.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: fxchain
	./fxchain

clean:
	rm -f fxchain
