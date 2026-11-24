CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -I../stream

all: spatial

spatial: src/spatial.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: spatial
	./spatial

clean:
	rm -f spatial
