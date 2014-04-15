FC = gfortran
FCFLAGS = -g -O3 -Wall

TARGETS = cgle.exe

all:$(TARGETS)

cgle.exe:cgle.f90
	$(FC) $(FCFLAGS) -o $@ $<

.PHONY:clean
clean:
	rm -r *.mod *.dat *.exe
