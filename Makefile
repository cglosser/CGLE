FC = gfortran
FFLAGS = -g -O3 -Wall -march=native

TARGETS = cgle.exe

all:$(TARGETS)

%.exe:%.f90
	$(FC) $(FFLAGS) -o $@ $<

.PHONY:clean
clean:
	rm -r *.mod *.dat *.exe
