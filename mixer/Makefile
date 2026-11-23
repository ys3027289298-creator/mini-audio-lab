CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

# NOGROUP=1 builds the pre-change unity-gain path (MIX_GAIN_GROUP=0). Its
# binary has a distinct name so the two configurations cannot go stale, and
# its test output is byte-identical to the original mixer.
ifeq ($(NOGROUP),1)
CFLAGS += -DMIX_GAIN_GROUP=0
BIN := mixer_nogroup
else
BIN := mixer
endif

SRC := src/mixer.c src/gain.c

all: $(BIN)

$(BIN): $(SRC) src/gain.h
	$(CC) $(CFLAGS) $(SRC) -o $@ -lm

test: $(BIN)
	./$(BIN)

ifeq ($(NOGROUP),1)
evidence:
	@echo "evidence unavailable: the gain structure is compiled out (NOGROUP=1)" >&2; exit 1
report:
	@echo "report unavailable: the gain structure is compiled out (NOGROUP=1)" >&2; exit 1
else
evidence: mixer
	./mixer --evidence

report: mixer
	@mkdir -p out
	./mixer --evidence | tee out/evidence.txt
endif

clean:
	rm -f mixer mixer_nogroup
	rm -rf out

.PHONY: all test evidence report clean
