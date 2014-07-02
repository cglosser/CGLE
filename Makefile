FC = gfortran
ENABLE_VISUALIZATION=true
#FFLAGS = -g -fcheck=all
FFLAGS =-O3 -Wall -march=native -ffast-math -cpp

EXECUTABLE = cgle.exe
TARGETS = latticeconst.o simparam.o cgle.o

ifdef ENABLE_VISUALIZATION
TARGETS += cglevis.o
FFLAGS += -D VISUALIZATION=1
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS = $(shell pkg-config --libs plplotd-f95) 
endif

COMPILE=$(FC) $(FFLAGS) -c

.PHONY:all
all:$(EXECUTABLE) Makefile

$(EXECUTABLE):$(TARGETS)
	$(FC) $(LIBS) -o $@ $^

cgle.o: cgle.f95 simparam.mod d2q9const.mod
	$(COMPILE) $<

cglevis.o cglevis.mod: cglevis.f95 simparam.mod
	$(COMPILE) $<

latticeconst.o d2q9const.mod d2q7const.mod: latticeconst.f95 simparam.mod
	$(COMPILE) $<

simparam.o simparam.mod: simparam.f95
	$(COMPILE) $<

.PHONY:clean
clean:
	rm -r *.mod *.dat *.exe *.o
