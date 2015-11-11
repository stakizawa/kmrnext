# Test programs in this directory are supported only on the following
# environment.
# - CentOS Linux 6.7
#   - gcc 4.4.7
#   - google test 1.7.0
# - OSX 10.11
#   - Apple LLVM version 7.0.0
#   - google test 1.7.0
#
include ../config.make

CXXFLAGS = -DDEBUG -g -O2 -Wall -std=c++98 -fopenmp
LIB      = ../src/libkmrnext.a

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
endif

TEST_OBJS = test_key.o test_data.o test_data_store.o
TEST_MAIN = test_all

# google test path
# The headers and libraries should be located under the directory as follows.
#   GTEST_HOME/
#     - include/
#       - gest/
#     - lib/
#       - libgtest.a
#       - libgtest_main.a  (not used in this test)
GTEST_HOME = $(HOME)/lib/googletest-1.7.0

.Phony: all
all: $(TEST_MAIN)

.cpp.o:
	$(CXX) $(DEFS) $(CXXFLAGS) -I$(GTEST_HOME)/include -I../src -c $<

test_all: test_main.o $(TEST_OBJS)
	$(CXX) -o $@ $^ $(LIB) $(GTEST_HOME)/lib/libgtest.a $(DEPLIB) -lpthread

.PHONY: clean
clean:
	-rm -rf $(TEST_MAIN) $(TESTS) *.o
