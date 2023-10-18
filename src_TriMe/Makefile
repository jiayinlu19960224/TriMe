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

# List of executables
execs= parallel_meshing_2d_main parallel_meshing_2d_main_fixed_point

# Makefile rules
all:
	$(MAKE) -C $(v2d)
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

parallel_meshing_2d_main: parallel_meshing_2d_main.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++

parallel_meshing_2d_main_fixed_point: parallel_meshing_2d_main_fixed_point.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++


%.o:%.cc
	$(cxx) $(cflags) $(iflags) -c $<


clean:
	rm -f $(execs)

.PHONY: all clean executables depend
