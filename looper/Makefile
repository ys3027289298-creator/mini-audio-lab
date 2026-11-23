CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: looper

looper: src/looper.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: looper
	./looper

plot: looper
	python3 -I tools/plot.py

clean:
	rm -f looper seam_wave.csv seam_diff.csv leakage_spec.csv seam_report.svg
