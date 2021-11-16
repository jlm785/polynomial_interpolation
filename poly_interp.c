/* poly_interp.f -- translated by f2c (version 20160102).
  from an old fortran version of poly_interp.f90.
  Cleaned by hand to not require libf2c
*/

// #include "f2c.h"

#include <stdlib.h>

/* >  Performs a Lagrange interpolation of a function and gives its */
/* >  Derivatives. */
/* >  Based on D. B. Hunter, The Computer Journal, 3, 270 (1961) */
/* >  with the small difference algorithm of Numerical Recipes Eq. 3.1.5 */
/* > */
/* >  \author       Jose Luis Martins, Carlos Loia Reis */
/* >  \version      6.0.6 */
/* >  \date         5 July 2021, 12 November 2021. */
/* >  \copyright    GNU Public License v3 */
/* Subroutine */ int poly_interp(double *y, double *dy, double *
	xin, double *yin, int *n, int *nd)
{
    /* System generated locals */
    int cmi_dim1, cmi_offset, dmi_dim1, dmi_offset, i1, i2, i3;
    double d1, d2;

    /* Local variables */
    static int i, k, m;
    static double dm, dp;
    static int ns, kmax;
    double *cmi=(double *) malloc(sizeof(double)*((*n+1) * (*nd+1)));
    double *dmi=(double *) malloc(sizeof(double)*((*n+1) * (*nd+1)));
    double *xnum=(double *) malloc(sizeof(double)*(*nd+1));

    /* Parameter adjustments */
    dmi_dim1 = *nd + 1;
    dmi_offset = 0;
    dmi -= dmi_offset;
    cmi_dim1 = *nd + 1;
    cmi_offset = 0;
    cmi -= cmi_offset;

    /* Function Body */
    ns = 0;
    i1 = *n;
    for (i = 0; i <= i1; ++i) {
	if ((d1 = xin[i], abs(d1)) < (d2 = xin[ns], abs(d2))) {
	    ns = i;
	}
    }
    y[0] = yin[ns];
    --ns;
    i1 = *n;
    for (i = 0; i <= i1; ++i) {
	cmi[i * cmi_dim1] = yin[i];
	dmi[i * dmi_dim1] = yin[i];
    }
    if (*nd > 0) {
	i1 = *nd;
	for (k = 1; k <= i1; ++k) {
	    y[k] = 0.;
	}
	i1 = *n;
	for (i = 0; i <= i1; ++i) {
	    i2 = *nd;
	    for (k = 1; k <= i2; ++k) {
		cmi[k + i * cmi_dim1] = 0.;
		dmi[k + i * dmi_dim1] = 0.;
	    }
	}
    }
    i1 = *n;
    for (m = 1; m <= i1; ++m) {
	kmax = m;
	if (kmax > *nd) {
	    kmax = *nd;
	}
	i2 = *n - m;
	for (i = 0; i <= i2; ++i) {
	    dm = xin[i];
	    dp = xin[i + m];
	    i3 = kmax;
	    for (k = 0; k <= i3; ++k) {
		xnum[k] = (dmi[k + i * dmi_dim1] - cmi[k + (i + 1) *
			cmi_dim1]) / (dp - dm);
	    }
	    dmi[i * dmi_dim1] = dp * xnum[0];
	    cmi[i * cmi_dim1] = dm * xnum[0];
	    if (*nd > 0) {
		i3 = kmax;
		for (k = 1; k <= i3; ++k) {
		    dmi[k + i * dmi_dim1] = dp * xnum[k] - k * xnum[k - 1];
		    cmi[k + i * cmi_dim1] = dm * xnum[k] - k * xnum[k - 1];
		}
	    }
	}
	if ((ns << 1) + 1 < *n - m) {
	    i2 = kmax;
	    for (k = 0; k <= i2; ++k) {
		dy[k] = cmi[k + (ns + 1) * cmi_dim1];
	    }
	} else {
	    i2 = kmax;
	    for (k = 0; k <= i2; ++k) {
		dy[k] = dmi[k + ns * dmi_dim1];
	    }
	    --ns;
	}
	i2 = kmax;
	for (k = 0; k <= i2; ++k) {
	    y[k] += dy[k];
	}
    }
    return 0;
} /* poly_interp__ */

