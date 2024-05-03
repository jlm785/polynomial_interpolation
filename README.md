# polynomial_interpolation
Lagrange polynomial interpolation with derivatives and error indication

Lagrange polynomial interpolation is sufficiently simple that for low fixed order of the polynomial it can be
coded essentially in a single line.

Here you can find a very short Fortran90 or C code that performs a Lagrange interpolation from an arbitrary set of data points
that includes the derivatives of the interpolated function and an indication of the quality of the interpolation.

This code is accompanied by an auxiliary code to interpolate from one grid to another (with error trapping and warnings),
testing examples plus a short description of the algorithm.

To run a simple test just type

$F90  test_poly.f90  poly_interp.f90
a.out

$F90  test_grid_interp.f90 grid_interp.f90 poly_interp.f90 sort.f90
a.out

where $F90 is any fortran compiler (ifort, gfortran,...).

OR

$CC test_poly.c  poly_interp.c
a.out

where $CC is any C compiler (gcc,...)
