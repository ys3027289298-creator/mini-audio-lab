CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -I../stream

all: mixer

mixer: src/mixer.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: mixer
	./mixer

clean:
	rm -f mixer
