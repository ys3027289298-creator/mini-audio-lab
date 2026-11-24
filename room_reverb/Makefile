PYTHON ?= python3
SRC := src

.PHONY: all run test clean

all: run

run:
	cd $(SRC) && $(PYTHON) experiment.py

# Quick physics sanity checks (image positions, reflection counts, Sabine limits).
test:
	cd $(SRC) && $(PYTHON) test_room.py

clean:
	rm -f artifacts/*.png artifacts/*.csv
	rm -rf src/__pycache__
