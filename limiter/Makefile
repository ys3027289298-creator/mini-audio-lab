CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: limiter

limiter: src/limiter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: limiter
	./limiter

clean:
	rm -f limiter
