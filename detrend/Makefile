CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: detrend bench

detrend: src/detrend.c src/detrend.h
	$(CC) $(CFLAGS) -DDETREND_SELFTEST src/detrend.c -o $@ -lm

bench: src/bench.c src/detrend.c src/detrend.h
	$(CC) $(CFLAGS) src/bench.c src/detrend.c -o $@ -lm

test: detrend
	./detrend

run: bench
	./bench

clean:
	rm -f detrend bench
