# Load the common configuration file 
include ../config.mk

v2d=../src_voro
tm2d=../src_TriMe
iflags=-I$(v2d) -I$(tm2d)
lflags=-L$(v2d) -L$(tm2d)

# Corresponding source code files
src=$(patsubst %.o,%.cc,$(objs))

# List of executables
execs= primitive_shape_meshing_sizing_function primitive_shape_meshing custom_poker_meshing parallel_meshing_2d_main parallel_meshing_2d_main_fixed_point parallel_meshing_2d_main_mapping

# Makefile rules
all:
	$(MAKE) -C $(v2d)
	$(MAKE) -C $(tm2d)
	$(MAKE) executables

executables: $(execs)

depend: 
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

-include Makefile.dep

primitive_shape_meshing_sizing_function: primitive_shape_meshing_sizing_function.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++ -ltrime++

primitive_shape_meshing: primitive_shape_meshing.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++ -ltrime++

custom_poker_meshing: custom_poker_meshing.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++ -ltrime++

parallel_meshing_2d_main: parallel_meshing_2d_main.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++ -ltrime++

parallel_meshing_2d_main_fixed_point: parallel_meshing_2d_main_fixed_point.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++ -ltrime++

parallel_meshing_2d_main_mapping: parallel_meshing_2d_main_mapping.cc $(objs)
	$(cxx) $(cflags) $(iflags) $(lflags) -o $@ $^ -lvoro++ -ltrime++

%.o:%.cc
	$(cxx) $(cflags) $(iflags) -c $<


clean:
	rm -f $(execs)

.PHONY: all clean executables depend
