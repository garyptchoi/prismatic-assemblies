# Include the configuration file
-include config.mk


# Lists of files to be built
heads=prismatic.hpp
objs=prismatic.o
src=$(patsubst %.o,%.cpp,$(objs) $(mg_objs))
execs=dofvscst

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm $(objs) $(execs)

%.o: %.cpp
	$(cxx) $(cflags) $(iflags) -c $< -I your-path-to-SuiteSparse/SuiteSparse-5.7.1/Include/ -L your-path-to-SuiteSparse/SuiteSparse-5.7.1/lib/ -lcholmod -lspqr

dofvscst: dofvscst.cpp $(objs)
				$(cxx) $(cflags) -o $@ $^ -I your-path-to-SuiteSparse/SuiteSparse-5.7.1/Include/ -L your-path-to-SuiteSparse/SuiteSparse-5.7.1/lib/ -lcholmod -lspqr

.PHONY: clean depend


# typing 'make' will invoke the first target entry in the file
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
#default: prismatic

# To create the executable file count we need the object files
# countwords.o, counter.o, and scanner.o:
#
#main:  main.cpp
#	$(CC) $(CFLAGS) -std=c++11 main.cpp
