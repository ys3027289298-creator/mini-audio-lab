CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: tfres

tfres: src/tfres.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: tfres
	./tfres

clean:
	rm -f tfres
