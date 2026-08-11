#include "crtbp.h"  // CRTBP2D class

#include "math_utils.h"  // astro::sqr

#include <cmath>  // std::sqrt

void CRTBP2D::fun(double t, const double *y, double *dydt, void *par) const
{
    (void)t;

    const auto  *p  = static_cast<const CRTBP2DParams *>(par);
    const double mu = p->mu;

    const double r1 = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
    const double r2 = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));

    const double r1_3 = 1.0 / (r1 * r1 * r1);
    const double r2_3 = 1.0 / (r2 * r2 * r2);

    dydt[0] = y[2];
    dydt[1] = y[3];

    dydt[2] = 2.0 * y[3] + y[0] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;

    dydt[3] = -2.0 * y[2] + y[1] * (1.0 - (1.0 - mu) * r1_3 - mu * r2_3);
}

void CRTBP2D::varfun(double t, const double *y, double *dydt, void *par) const
{
    (void)t;

    const auto  *p  = static_cast<const CRTBP2DParams *>(par);
    const double mu = p->mu;

    const double r1 = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
    const double r2 = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));

    const double r1_3 = 1.0 / (r1 * r1 * r1);
    const double r2_3 = 1.0 / (r2 * r2 * r2);
    const double r1_5 = r1_3 / astro::sqr(r1);
    const double r2_5 = r2_3 / astro::sqr(r2);

    // Equations of motion.
    dydt[0] = y[2];
    dydt[1] = y[3];

    dydt[2] = 2.0 * y[3] + y[0] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;

    dydt[3] = -2.0 * y[2] + y[1] * (1.0 - (1.0 - mu) * r1_3 - mu * r2_3);

    // Second derivatives of the effective potential.
    const double O_xx = 1.0 - (1.0 - mu) * r1_3 - mu * r2_3 + 3.0 * (1.0 - mu) * astro::sqr(y[0] + mu) * r1_5 +
                        3.0 * mu * astro::sqr(y[0] - 1.0 + mu) * r2_5;

    const double O_xy = 3.0 * (1.0 - mu) * (y[0] + mu) * y[1] * r1_5 + 3.0 * mu * (y[0] - 1.0 + mu) * y[1] * r2_5;

    const double O_yy = 1.0 - (1.0 - mu) * r1_3 - mu * r2_3 + 3.0 * (1.0 - mu) * astro::sqr(y[1]) * r1_5 +
                        3.0 * mu * astro::sqr(y[1]) * r2_5;

    // Variational equations.
    dydt[4] = y[6];
    dydt[5] = y[7];

    dydt[6] = O_xx * y[4] + O_xy * y[5] + 2.0 * y[7];
    dydt[7] = O_xy * y[4] + O_yy * y[5] - 2.0 * y[6];
}