CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: downmix

downmix: src/downmix.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: downmix
	./downmix
	cmp bypass_in.f32 bypass_out.f32
	@echo "bypass_cmp=identical"

clean:
	rm -f downmix bypass_in.f32 bypass_out.f32 \
	      downmix_in.f32 downmix_lt.f32 downmix_rt.f32
