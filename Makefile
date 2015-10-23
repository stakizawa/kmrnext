CC = g++
CFLAGS = -Wall -O2

.Phony: all
all: sample

sample: sample.cpp
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean
clean:
	-rm -rf sample

