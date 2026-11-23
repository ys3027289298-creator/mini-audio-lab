CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: doppler

doppler: src/doppler.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: doppler
	./doppler

# Rebuilds the CSV/JSON under ../artifacts/doppler from the same source.
report: tools/report.c src/doppler.c
	$(CC) $(CFLAGS) -Wno-unused-function tools/report.c -o /tmp/doppler_report -lm
	mkdir -p ../artifacts/doppler
	/tmp/doppler_report ../artifacts/doppler

clean:
	rm -f doppler
