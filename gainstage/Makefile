CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: gainstage

gainstage: src/gainstage.c ../limiter/src/limiter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: gainstage
	./gainstage

clean:
	rm -f gainstage
