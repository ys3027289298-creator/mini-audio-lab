CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: meter lab

meter: src/meter.c src/meter.h
	$(CC) $(CFLAGS) src/meter.c -o $@ -lm

lab: src/lab.c src/meter.h
	$(CC) $(CFLAGS) src/lab.c -o $@ -lm

report: lab
	./lab
	./lab --no-unify

test: meter lab
	./meter
	./lab
	./lab --no-unify
	cmp out/export_1024.raw out/playback_1024.raw
	cmp out/export_1000.raw out/playback_1000.raw
	cmp out/nounify/export_1024.raw out/nounify/playback_1024.raw
	cmp out/nounify/export_1000.raw out/nounify/playback_1000.raw

clean:
	rm -f meter lab
	rm -rf out

.PHONY: all report test clean
