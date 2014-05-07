FC = gfortran
FFLAGS = -g -O3 -Wall -march=native
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS = $(shell pkg-config --libs plplotd-f95)

EXECUTABLE = cgle.exe
TARGETS = latticeconst.o simparam.o cglevis.o cgle.o

all:$(EXECUTABLE)

$(EXECUTABLE):$(TARGETS)
	$(FC) $(LIBS) -o $@ $^

%.o:%.f95
	$(FC) $(FFLAGS) -c $<

.PHONY:clean
clean:
	rm -r *.mod *.dat *.exe *.o
