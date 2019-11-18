FC=h5pfc
CC=gcc
FCFLAGS=-O3 -Wall
CCFLAGS=-O3 -Wall 
LDFLAGS=

EXE=heat_mpi
OBJS=main.o heat_mod.o core.o setup.o utilities.o io.o

all: $(EXE)

core.o: core.F90 heat_mod.o
utilities.o: utilities.F90 heat_mod.o
io.o: io.F90 heat_mod.o 
setup.o: setup.F90 heat_mod.o utilities.o io.o
main.o: main.F90 heat_mod.o core.o io.o setup.o utilities.o

$(EXE): $(OBJS) 
	$(FC) $(FCFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LIBS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CCFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-/bin/rm -f $(EXE) a.out *.o *.mod *.png *~
