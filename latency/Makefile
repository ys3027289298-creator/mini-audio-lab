CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: latency

latency: src/latency.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: latency
	./latency

clean:
	rm -f latency
