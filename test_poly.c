/* test_poly.f -- translated by f2c (version 20160102).
   from an old fortran version of test_poly.f90
   Cleaned to not require libf2c
*/
#include <stdio.h>

/* Table of constant values */

static int c__0 = 0;
static int c__1 = 1;

/* Main program */ int main(int argc, char **argv)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Local variables */
    static int j, k, n;
    static double x, y[11];
    static int nd;
    static double dy[11], xin[21], yin[21], yex[11];
    extern /* Subroutine */ int poly_interp(double *, double *,
	    double *, double *, int *, int *), fexact(
	    double *, double *, int *);


    n = 5;
    nd = n;
    i1 = n;
    for (j = 0; j <= i1; ++j) {
	xin[j] = j + 10.f;
	fexact(&xin[j], yex, &c__0);
	yin[j] = yex[0];
    }
    x = 12.3f;
    fexact(&x, yex, &nd);
    i1 = n;
    for (j = 0; j <= i1; ++j) {
	xin[j] -= x;
    }
    poly_interp(y, dy, xin, yin, &n, &nd);
    i1 = nd;
    printf("%5i \n",i1);
    for (k = 0; k <= i1; ++k) {
     printf("%14.8lf %14.8lf %14.8lf %14.8lf \n",y[k],yex[k],y[k]-yex[k],dy[k]);
    }
    return 0;
} /* MAIN__ */


/* fexact_sqrt.f -- translated by f2c (version 20160102).
   from an old fortran version of fexact.f90
   Cleaned to not require libf2c
*/

#include <math.h>

/* Subroutine */ int fexact(double *x, double *y, int *nd)
{
    /* System generated locals */
    int i1;
    double d1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int k;
    static double factor;

    y[0] = sqrt(*x);
    factor = 1.;
    if (*nd > 0) {
	i1 = *nd;
	for (k = 1; k <= i1; ++k) {
	    factor = factor * (3. - (k << 1)) / 2.;
	    d1 = .5 - k * 1.;
	    y[k] = factor * pow(*x, d1);
	}
    }
    return 0;
} /* fexact_ */


