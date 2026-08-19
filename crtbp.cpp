#include "crtbp.h"  // CRTBP2D class

#include "math_utils.h"  // astro::sqr

#include <cmath>    // std::sqrt
#include <iomanip>  // std::scientific, std::setprecision, std::setw
#include <ostream>  // std::ostream

CRTBP2D::CRTBP2D(double mu)
{
    setName("Planar CRTBP");
    setNVar(4);

    param_.mu = mu;
}

void CRTBP2D::InertialToCRTBP(const astro::State &state, double a2, double n)
{
    const double mu = param_.mu;
    double      *y  = getY();

    const double xi     = state.r.x;
    const double eta    = state.r.y;
    const double xiDot  = state.v.x;
    const double etaDot = state.v.y;

    // Dimensionless barycentric rotating position.
    y[0] = xi / a2 - mu;
    y[1] = eta / a2;

    // Dimensionless barycentric rotating velocity.
    y[2] = xiDot / (n * a2) + eta / a2;
    y[3] = etaDot / (n * a2) - xi / a2;
}

void CRTBP2D::GetInitialCondition(double a, double e)
{
    const double mu = param_.mu;
    double      *y  = getY();

    // Initial position at periapsis.
    y[0] = a * (1.0 - e) - mu;
    y[1] = 0.0;

    // Initial velocity perpendicular to the x-axis.
    y[2] = 0.0;
    y[3] = std::sqrt((1.0 - mu) / a * (1.0 + e) / (1.0 - e)) - a * (1.0 - e);
}

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

void CRTBP2D::printState(std::ostream &os, double t, const double *y) const
{
    constexpr int W         = 18;
    constexpr int precision = 10;

    os << std::scientific << std::showpos << std::setprecision(precision);

    os << std::setw(W) << t;
    for (std::size_t i = 0; i < getNVar(); ++i) {
        os << std::setw(W) << y[i];
    }

    os << '\n';
}