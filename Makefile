CC = g++
CFLAGS = -Wall -O2

.PHONY: all clean

all: sample

sample: sample.cpp
	$(CC) $(CFLAGS) -o $@ $^

clean:
	-rm -rf sample

