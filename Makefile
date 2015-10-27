CC = g++
CFLAGS = -Wall -O2

HEADERS = kmrnext.hpp
OBJS    = kmrnext.o
TESTS   = test00-basic test01-mng-ds

.Phony: all
all: $(TESTS)

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) -c $<

test00-basic: $(OBJS) test00-basic.o
	$(CC) -o $@ $^

test01-mng-ds: $(OBJS) test01-mng-ds.o
	$(CC) -o $@ $^

.PHONY: clean
clean:
	-rm -rf $(OBJS)
	-rm -rf $(TESTS) *.o

