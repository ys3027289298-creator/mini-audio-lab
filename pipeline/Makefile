CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -I../stream -I../resampler/src -I../mixer/src -I../spatial/src

ALL = pipeline_legacy pipeline_unified

all: $(ALL)

pipeline_legacy: pipeline.c
	$(CC) $(CFLAGS) -DMA_UNIFY=0 $< -o $@ -lm

pipeline_unified: pipeline.c
	$(CC) $(CFLAGS) -DMA_UNIFY=1 $< -o $@ -lm

test: $(ALL)
	@mkdir -p /tmp/ma_demo_legacy /tmp/ma_demo_unified
	./pipeline_legacy demo /tmp/ma_demo_legacy
	./pipeline_unified demo /tmp/ma_demo_unified

clean:
	rm -f $(ALL)
