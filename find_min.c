#include "include/fields_nate.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

double fn1 (double z_in, void * params) {
    (void)(params); /* avoid unused parameter warning */
    double x, y, z, fx, fy, fz, u, t, freq;
    x = 0.0;
    y = 0.0;
    z = z_in;
    t = 0.0;
    freq = 0.0;
    force_(&x, &y, &z, &fx, &fy, &fz, &u, &t, &freq);
    return u;
}

int main (void) {
    int status;
    int iter = 0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;

    double min = -1.45;
    double lowb = -1.49999;
    double upb = -1.3;

    gsl_function F;

    F.function = &fn1;
    F.params = 0;

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    gsl_min_fminimizer_set (s, &F, min, lowb, upb);

    printf ("using %s method\n",
          gsl_min_fminimizer_name (s));

    printf ("%5s [%9s, %9s] %9s %9s\n",
          "iter", "lower", "upper", "min", "err(est)");

    printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
          iter, lowb, upb,
          min, upb - lowb);

    do {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        min = gsl_min_fminimizer_x_minimum (s);
        lowb = gsl_min_fminimizer_x_lower (s);
        upb = gsl_min_fminimizer_x_upper (s);

        status = gsl_min_test_interval (lowb, upb, 0.0, 1e-6);

        if (status == GSL_SUCCESS) {
            printf ("Converged:\n");
        }

        printf ("%5d [%.7f, %.7f] %.15e %.15e\n", iter, lowb, upb, min, upb - lowb);
        printf("%.15e\n", gsl_min_fminimizer_f_minimum(s));
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free (s);

    return status;
}