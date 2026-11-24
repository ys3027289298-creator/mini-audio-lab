CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: limiter

limiter: src/limiter.c ../common/sr.h
	$(CC) $(CFLAGS) -I../common $< -o $@ -lm

test: limiter
	./limiter

clean:
	rm -f limiter
