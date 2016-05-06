all:bin/program_train_emission bin/program_external_information bin/program_HMM_modelling


LINK.o = $(LINK.cpp)
CXXFLAGS = -Wall -O2 -std=gnu++0x -fopenmp
CXX = g++
objs0 = program_train_emission
objs1 = program_external_information 
objs2 = program_HMM_modelling
objs0 := $(addsuffix .o, $(objs0))
objs0 := $(addprefix src/, $(objs0))
objs1 := $(addsuffix .o, $(objs1))
objs1 := $(addprefix src/, $(objs1))
objs2 := $(addsuffix .o, $(objs2))
objs2 := $(addprefix src/, $(objs2))
CPPFLAGS = -I$(HOME)/boost/include -I.

bin/program_train_emission : $(objs0)
	$(CXX) $(objs0) -fopenmp -o $@ 

bin/program_external_information: $(objs1)
	$(CXX) $(objs1) -fopenmp -o $@ 

bin/program_HMM_modelling: $(objs2)
	$(CXX) $(objs2) -fopenmp -o $@ 



clean:
	$(RM) bin/program_train_emission bin/program_external_information bin/program_HMM_modelling $(objs)
	
