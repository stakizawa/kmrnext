CC = g++
CFLAGS = -Wall -O2

HEADERS = kmrnext.hpp
OBJS = kmrnext.o

.Phony: all
all: sample

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) -c $<

sample: $(OBJS) sample.o
	$(CC) -o $@ $^

.PHONY: clean
clean:
	-rm -rf $(OBJS)
	-rm -rf sample *.o

