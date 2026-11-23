CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

# NOLOCATE=1 builds the original occurrence-only detector (no clip intervals).
# Its binary has a distinct name so the two configurations cannot go stale.
ifeq ($(NOLOCATE),1)
CFLAGS += -DMIX_CLIP_LOCATE=0
BIN := mixer_nolocate
else
BIN := mixer
endif

all: $(BIN)

$(BIN): src/mixer.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: $(BIN)
	./$(BIN)

ifeq ($(NOLOCATE),1)
evidence:
	@echo "evidence unavailable: rebuild without NOLOCATE=1" >&2; exit 1
else
evidence: mixer
	./mixer --evidence
endif

clean:
	rm -f mixer mixer_nolocate
