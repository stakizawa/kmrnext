include ../config.make

LIB    = ../src/libkmrnext.a
PROGS  = map-basic map-scale\
	 nicam-letkf-pattern-cell nicam-letkf-pattern-region \
	 nicam-letkf-pattern-mpi
DEPLIB = $(KMRDIR)/lib/libkmr.a

.Phony: all
all: $(PROGS)

map-basic: map-basic.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

map-scale: map-scale.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

nicam-letkf-pattern-cell: nicam-letkf-pattern-cell.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

nicam-letkf-pattern-region: nicam-letkf-pattern-region.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

nicam-letkf-pattern-mpi: nicam-letkf-pattern-mpi.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) $(DEPLIB)

.cpp.o:
ifneq ($(TARGET), KMR)
	$(error KMR backend should be used)
endif
	$(CXX) $(DEFS) $(CXXFLAGS) -I../src -I$(KMRDIR)/include -c $<

.PHONY: clean
clean:
	-rm -rf $(PROGS) *.o *.d

.PHONY: distclean
distclean: clean
