CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: specquant

specquant: src/specquant.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: specquant
	./specquant

clean:
	rm -f specquant
