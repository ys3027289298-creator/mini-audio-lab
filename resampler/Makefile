CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: resampler

resampler: src/resampler.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: resampler
	./resampler

measure: aa_measure
	@mkdir -p doc
	./aa_measure | tee doc/measurements.md
	python3 tools/plot_freq.py doc/freq_response.csv doc/freq_response.svg

aa_measure: tools/measure.c src/resampler.c
	$(CC) $(CFLAGS) tools/measure.c -o $@ -lm

clean:
	rm -f resampler aa_measure

.PHONY: all test measure clean
