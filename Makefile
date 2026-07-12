PROG = CageBreak	

SRCS = 	cagebreak.f checkneigh.f io1.f io2.f neighboursAB.f timers.f

OBJS = 	cagebreak.o checkneigh.o io1.o io2.o neighboursAB.o timers.o

LIBS =	

FC = gfortran
F90 = gfortran
F90FLAGS=  -v -ffixed-line-length-132 -fbounds-check -O3 -march=native -ftree-vectorize -mcmodel=large -ftree-vectorize -fopt-info-vec -fipa-cp -finline-functions
FFLAGS=  -ffixed-line-length-132 -fbounds-check  -O3 -march=native -ftree-vectorize  -mcmodel=large  -fopt-info-vec -fipa-cp -finline-functions 
LDFLAGS = $(F90FLAGS)
#FC =ifort
#F90FLAGS= -132 -heap-arrays -Vaxlib -g -CB 
#FFLAGS= -132 -heap-arrays -Vaxlib -g -CB 
#LDFLAGS = $(F90FLAGS)

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *.oo

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

cagebreak.o: system.h counter.h criteria.h pos.h switch.h
checkneigh.o: system.h counter.h pos.h switch.h
io1.o: counter.h switch.h 
io2.o: system.h pos.h counter.h 
neighboursAB.o: system.h criteria.h counter.h switch.h pos.h 
timers.o: system.h switch.h

