#FFLAGS = -O3 -fbounds-check -ffree-form -ffree-line-length-none -fno-backslash -freal-4-real-8 
FFLAGS = -O3 -freal-4-real-8 -fPIC
# For quad precision:
#FFLAGS = -O3 -fbounds-check -ffree-form -ffree-line-length-none -fno-backslash -freal-4-real-16
# For debug:
#FFLAGS = -g -fbounds-check -ffree-form -ffree-line-length-none -fno-backslash -freal-4-real-8 

LAPACK_LIB = /usr/lib

shared: genquad.f90 quadpack.f90
	gfortran -shared $(FFLAGS) genquad.f90 -o genquad.so
	gfortran -shared $(FFLAGS) quadpack.f90 -o quadpack.so

genquad.o : genquad.f90 quadpack.o
	gfortran -c $(FFLAGS) genquad.f90

test : test.o genquad.o
	gfortran $(FFLAGS) -o test test.o genquad.o quadpack.o -L$LAPACK_LIB -llapack
	./test
	rm test

test.o : test.f90 genquad.o
	gfortran -c $(FFLAGS) test.f90

quadpack.o : quadpack.f90
	gfortran -c $(FFLAGS) quadpack.f90

clean:
	rm -f *.o *.mod test *.plt *.out *.dat
