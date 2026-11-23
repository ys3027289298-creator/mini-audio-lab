CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
# MIXER_SMOOTH_RAMP=0 restores the legacy linear ramp bit-for-bit.
RAMP_DEF ?= -DMIXER_SMOOTH_RAMP=1

all: mixer

mixer: src/mixer.c
	$(CC) $(CFLAGS) $(RAMP_DEF) $< -o $@ -lm

test: mixer
	./mixer

clean:
	rm -f mixer
