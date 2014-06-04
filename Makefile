FC = gfortran
FFLAGS = -g -O3 -Wall -march=native -ffast-math
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS = $(shell pkg-config --libs plplotd-f95)

EXECUTABLE = cgle.exe
TARGETS = latticeconst.o simparam.o cglevis.o cgle.o

COMPILE=$(FC) $(FFLAGS) -c

.PHONY:all
all:$(EXECUTABLE)

$(EXECUTABLE):$(TARGETS)
	$(FC) $(LIBS) -o $@ $^

cgle.o: cgle.f95 simparam.mod d2q9const.mod cglevis.mod
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
