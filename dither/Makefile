CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -pedantic
CPPFLAGS ?= -Iinclude

CORE = src/dither_core.c include/dither.h

all: dither dither_export

dither: src/dither_test.c $(CORE)
	$(CC) $(CPPFLAGS) $(CFLAGS) src/dither_test.c src/dither_core.c -o $@ -lm

dither_export: src/dither_export.c $(CORE)
	$(CC) $(CPPFLAGS) $(CFLAGS) src/dither_export.c src/dither_core.c -o $@ -lm

test: dither
	./dither

clean:
	rm -f dither dither_export
