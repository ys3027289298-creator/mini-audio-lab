CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: limiter

limiter: src/limiter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: limiter
	./limiter

compare: src/compare.c tools/plot_compare.py
	$(CC) $(CFLAGS) src/compare.c -o compare -lm
	mkdir -p out
	./compare
	python3 tools/plot_compare.py

clean:
	rm -f limiter compare
	rm -rf out
