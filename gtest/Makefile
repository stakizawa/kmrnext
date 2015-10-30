# Test programs in this directory is only supported on the following
# environment.
# - CentOS Linux 6.7
# - gcc 4.4.7
# - google test 1.7.0
#
include ../config.make

CXXFLAGS = -DDEBUG -g -O2 -Wall -std=c++98 -fopenmp
HEADERS  = kmrnext.hpp
LIB      = ../src/libkmrnext.a
TESTS    = test_key test_data

# google test path
GTEST_HOME = $(HOME)/lib/googletest-1.7.0

.Phony: all
all: $(TESTS)

.cpp.o:
	$(CXX) $(DEFS) $(CXXFLAGS) -I$(GTEST_HOME)/include -I../src -c $<

test_key: test_key.o
	$(CXX) -o $@ $^ $(LIB) $(GTEST_HOME)/lib/libgtest.a -lpthread

test_data: test_data.o
	$(CXX) -o $@ $^ $(LIB) $(GTEST_HOME)/lib/libgtest.a -lpthread

run_tests: $(TESTS)
	./test_key

.PHONY: clean
clean:
	-rm -rf $(TESTS) *.o *.d
