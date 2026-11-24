CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: latcal

latcal: src/latcal.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: latcal
	./latcal

clean:
	rm -f latcal
