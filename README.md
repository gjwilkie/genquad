# genquad
A routine that generates an arbitrarily-weighted quadrature scheme for Gaussian numerical integration. Written in Fortran 90 and requires the [LAPACK](www.netlib.org/lapack/) linear algebra library. 

You wish to approximate the integral:

![alt text](https://raw.githubusercontent.com/gjwilkie/markdownstuff/master/quadrature.gif)

where *x<sub>i</sub>* and *w<sub>i</sub>* are the grid points (a.k.a. "abscissae") and weights of the *N*-point quadrature scheme, respectively. Knowledge of the weighting function *g*(*x*) is contained in these arrays. This package generates these abscissae and weights for any given (relatively well-behaved) weighting function *g*.

An important feature of Gaussian quadrature is *spectral accuracy*. Instead of the error decreasing as a power of *N*, it decreases **exponentially* with *N*.

## Using genquad

1. Copy quadpack.f90 and genquad.f90 to your source directory and include the `genquad` module
  
  ```
  use genquad, only: get_quadrature_rule
  ```
2. When building your application, include your LAPACK library when linking. Also, include following in your Makefile:
  ```
  FFLAGS = -freal-r-real-8
  
  genquad.o : genquad.f90 quadpack.o
     gfortran -c $(FFLAGS) genquad.f90

  quadpack.o : quadpack.f90
     gfortran -c $(FFLAGS) quadpack.f90
  ```
3. Make sure the weighting function is defined as `external` in the calling routine.

## Testing genquad

To ensure integrals are acurately estimated and specral accuracy is maintained, you may wish to run the test suite. To do so:

1. Update Makefile with the location of your LAPACK library and, if necessary, your fortran compiler.
2. Run `make test`

## Reference

See Appendix C of [thesis](https://drive.google.com/open?id=0B5fJ4SuNBdN3TktSZTltYldKSDA) for more details.

