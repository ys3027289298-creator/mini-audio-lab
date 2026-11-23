CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: drift

drift: src/drift.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: drift
	./drift

clean:
	rm -f drift
