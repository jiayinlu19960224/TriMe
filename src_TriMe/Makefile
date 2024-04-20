# Load the common configuration file 
include ../config.mk

v2d=../src_voro
iflags=-I$(v2d)
lflags=-L$(v2d)

# List of objects to be built
objs=basic_calculation.o parallel_meshing_2d.o shape_2d.o custom_shape_2d.o \
	 adf_2d.o sizing_2d.o mesh_alg_2d.o mesh_alg_2d_dm.o \
	mesh_alg_2d_cvd.o \
	mesh_alg_2d_hybrid.o

# Corresponding source code files
src=$(patsubst %.o,%.cc,$(objs))

# Makefile rules
all: libtrime++.a 
	
depend: 
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

libtrime++.a: $(objs)
	rm -f libtrime++.a
	ar rs libtrime++.a $^

%.o:%.cc
	$(cxx) $(cflags) $(iflags) -c $<


help: Doxyfile $(SOURCE)
	doxygen Doxyfile

clean:
	rm -f $(objs) libtrime++.a

.PHONY: all help execs depend
