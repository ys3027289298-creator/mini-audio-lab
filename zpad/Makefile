CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: zpad

zpad: src/zpad.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: zpad
	./zpad

clean:
	rm -f zpad
