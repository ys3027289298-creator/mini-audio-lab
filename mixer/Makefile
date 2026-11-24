CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: mixer

mixer: src/mixer.c ../common/sr.h
	$(CC) $(CFLAGS) -I../common $< -o $@ -lm

test: mixer
	./mixer

clean:
	rm -f mixer
