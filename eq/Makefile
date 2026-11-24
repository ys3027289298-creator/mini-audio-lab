CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: eq

eq: src/eq.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: eq
	./eq

clean:
	rm -f eq
