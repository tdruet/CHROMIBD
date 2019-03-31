f90 = gfortran
prog=chromibd_1.2

a.out:  $(prog).o pedmosaic.o unrelated.o
	$(f90) $(prog).o pedmosaic.o unrelated.o
	mv a.out $(prog)

$(prog).o:      $(prog).f90 pedmosaic.o unrelated.o
	$(f90)  -c $(prog).f90

pedmosaic.o: pedmosaic.f90
	$(f90) -c pedmosaic.f90
unrelated.o: unrelated.f90
	$(f90) -c unrelated.f90
clean:
	rm *.o *.mod

