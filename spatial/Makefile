CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: spatial early_reflection

spatial: src/spatial.c
	$(CC) $(CFLAGS) $< -o $@ -lm

early_reflection: src/early_reflection.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: spatial early_reflection
	./spatial
	./early_reflection

clean:
	rm -f spatial early_reflection
