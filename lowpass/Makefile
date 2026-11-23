CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
PYTHON ?= python3

all: lowpass

lowpass: src/lowpass.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lowpass
	./lowpass

verify: test
	$(PYTHON) tools/verify.py

plots: test
	$(PYTHON) tools/plot_results.py

clean:
	rm -f lowpass sweep.csv cascade.csv pareto_front.png cascade.png

.PHONY: all test verify plots clean
