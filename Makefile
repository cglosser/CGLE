FC = gfortran
FFLAGS = -g -O3 -Wall -march=native
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS = $(shell pkg-config --libs plplotd-f95)

EXECUTABLE = cgle.exe
TARGETS = latticeconst.o simparam.o cglevis.o cgle.o

COMPILE=$(FC) $(FFLAGS) -c

all:$(EXECUTABLE)

$(EXECUTABLE):$(TARGETS)
	$(FC) $(LIBS) -o $@ $^

cgle.o: cgle.f95 cglevis.o simparam.o latticeconst.o
	$(COMPILE) $<

cglevis.o: cglevis.f95 simparam.o
	$(COMPILE) $<

latticeconst.o: latticeconst.f95 simparam.o
	$(COMPILE) $<

simparam.o: simparam.f95
	$(COMPILE) $<

#%.o:%.f95
	#$(FC) $(FFLAGS) -c $<

.PHONY:clean
clean:
	rm -r *.mod *.dat *.exe *.o
