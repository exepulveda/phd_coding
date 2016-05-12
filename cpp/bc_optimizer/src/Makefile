SRCS = Population.h \
		Random.cpp Population.cpp Crossover.cpp \
        Individual.cpp IntArrayIndividual.cpp \
        hdf5_utils.cpp \
        BlockModel.cpp DrawPoint.cpp BlockCavingProblem.cpp 

OBJS = $(SRCS:.cpp=.o)

ARMADILLO_HOME = /usr/local
HDF5_HOME = /usr/local
BOOST_HOME = /usr/local

INCLUDES = -I$(ARMADILLO_HOME)/include  -I$(HDF5_HOME)/include -I$(HDF5_HOME)/include

LIBRARIES = -L$(ARMADILLO_HOME)/lib  -L$(HDF5_HOME)/lib -L$(BOOST_HOME)/lib 

CC = g++
H5CC = h5c++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG) -std=c++11 -fopenmp 
LFLAGS = -Wall $(DEBUG) -std=c++11 -fopenmp $(LIBRARIES) -lz -lm -larmadillo -fopenmp -lboost_system -lboost_filesystem 

all: bc_mo_plain bc_so_plain

bc_mo_plain : $(OBJS) bc_mo_plain.o
	$(H5CC) $(LFLAGS) bc_mo_plain.o $(OBJS) -o bc_mo_plain 

bc_so_plain : $(OBJS) bc_so_plain.o
	$(H5CC) $(LFLAGS) bc_so_plain.o $(OBJS) -o bc_so_plain 

bc_mo_plain.o: bc_mo_plain.cpp $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

bc_so_plain.o: bc_so_plain.cpp $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

.cpp.o: $(SRCS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	\rm *.o bc_mo_plain bc_so_plain
