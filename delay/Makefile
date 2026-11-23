CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: delay

delay: src/delay.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: delay
	./delay

# measurements for the analysis figures
csv: delay
	mkdir -p ../artifacts/delay
	./delay > ../artifacts/delay/measurements.txt

clean:
	rm -f delay
