CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: peakpick

peakpick: src/peakpick.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: peakpick
	./peakpick

clean:
	rm -f peakpick
