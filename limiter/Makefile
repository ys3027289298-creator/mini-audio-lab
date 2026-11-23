CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: limiter

limiter: src/limiter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: limiter
	./limiter

# lookahead transient study: sweep + waveforms into ../artifacts/lookahead
study:
	$(CC) $(CFLAGS) exp/lookahead_study.c -o /tmp/lookahead_study -lm
	mkdir -p ../artifacts/lookahead
	/tmp/lookahead_study ../artifacts/lookahead
	python3 exp/plot_lookahead.py ../artifacts/lookahead

clean:
	rm -f limiter
