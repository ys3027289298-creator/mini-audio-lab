CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hop

hop: src/hop.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: hop
	./hop

clean:
	rm -f hop
