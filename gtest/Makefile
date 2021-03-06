# Test programs in this directory are supported only on the following
# environment.
# - CentOS Linux 6.7
#   - gcc 4.4.7
#   - google test 1.7.0
# - K computer
#   - fujitsu compiler (latest)
#   - google test 1.7.0
#
include ../config.make

LIB      = ../src/libkmrnext.a
TEST_OBJS = test_kmrnext.o test_key.o test_data.o test_data_store.o \
            test_data_store_view.o
TEST_MAIN = test_all

ifeq ($(SYSENV), GNULINUX)
    CXXFLAGS += --std=gnu++0x
endif

ifeq ($(SYSENV), K_FX)
    CXXFLAGS += -DGTEST_HAS_TR1_TUPLE=0
endif

ifeq ($(TARGET), KMR)
    DEPLIB += $(KMRDIR)/lib/libkmr.a
    TEST_OBJS += test_data_store_kmr.o
endif

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
ifeq ($(TARGET), KMR)
	$(CXX) $(DEFS) $(CXXFLAGS) -I$(GTEST_HOME)/include \
	-I$(KMRDIR)/include -I../src -c $<
else
	$(CXX) $(DEFS) $(CXXFLAGS) -I$(GTEST_HOME)/include -I../src -c $<
endif

test_all: test_main.o $(TEST_OBJS)
	$(CXX)  $(LDFLAGS) -o $@ $^ $(LIB) $(GTEST_HOME)/lib/libgtest.a \
	$(DEPLIB) -lpthread

run_tests: $(TEST_MAIN)
ifeq ($(TARGET), KMR)
	mpiexec --tag-output -n 2 ./$(TEST_MAIN)
else
	./$(TEST_MAIN)
endif

.PHONY: clean
clean:
	-rm -rf $(TEST_MAIN) $(TESTS) *.o *.d

-include *.d
