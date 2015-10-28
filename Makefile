CC = g++
CFLAGS = -Wall -MD -g

HEADERS = kmrnext.hpp
OBJS    = kmrnext.o
TESTS   = test00-basic test01-mng-ds

.Phony: all
all: $(TESTS)

test00-basic: $(OBJS) test00-basic.o
	$(CC) -o $@ $^

test01-mng-ds: $(OBJS) test01-mng-ds.o
	$(CC) -o $@ $^

.cpp.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	-rm -rf $(OBJS)
	-rm -rf $(TESTS) *.o *.d

-include *.d

