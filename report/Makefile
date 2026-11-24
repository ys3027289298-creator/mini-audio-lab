CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
COMMON = -I../common

all: sr_report

# Compile the real module sources with static linkage removed and
# conflicting symbols renamed, so the report links the actual DSP code.
rs.o: ../resampler/src/resampler.c ../common/sr.h
	$(CC) $(CFLAGS) $(COMMON) -Dstatic= -Dmain=rs_selftest -Dcheck=rs_check -Dfailures=rs_failures -Dpeak=rs_peak -c $< -o $@

mx.o: ../mixer/src/mixer.c ../common/sr.h
	$(CC) $(CFLAGS) $(COMMON) -Dstatic= -Dmain=mx_selftest -Dcheck=mx_check -Dfailures=mx_failures -c $< -o $@

lm.o: ../limiter/src/limiter.c ../common/sr.h
	$(CC) $(CFLAGS) $(COMMON) -Dstatic= -Dmain=lm_selftest -Dcheck=lm_check -Dfailures=lm_failures -Dcoef=lm_coef -c $< -o $@

sr_report: sr_report.c rs.o mx.o lm.o ../common/sr.h
	$(CC) $(CFLAGS) $(COMMON) sr_report.c rs.o mx.o lm.o -o $@ -lm

test: sr_report
	./sr_report

# Byte-identity proof: rebuild every module's self-test and diff its output
# against the golden files captured before the unification.
check-compat:
	@for m in resampler mixer limiter fft spatial; do \
		$(MAKE) -s -C ../$$m >/dev/null; \
		../$$m/$$m > /tmp/mal_now_$$m.out; \
		if cmp -s golden/$$m.out /tmp/mal_now_$$m.out; then \
			echo "$$m: byte-identical"; \
		else \
			echo "$$m: OUTPUT CHANGED"; diff golden/$$m.out /tmp/mal_now_$$m.out; exit 1; \
		fi; \
	done

clean:
	rm -f sr_report rs.o mx.o lm.o

.PHONY: all test check-compat clean
