CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

SRC = src/spatial.c src/atten.c
HDR = src/atten.h

all: spatial

spatial: $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(SRC) -o $@ -lm

report: spatial
	./spatial --report

test: spatial
	./spatial
	./spatial --report

clean:
	rm -f spatial
	rm -rf out

.PHONY: all report test clean
