#include "crtbp.h"       // CRTBP2D class
#include "io.h"          // Common numerical output formatting
#include "math_utils.h"  // astro::sqr

#include <cmath>      // std::sqrt
#include <iomanip>    // std::scientific, std::setprecision, std::setw
#include <ostream>    // std::ostream
#include <stdexcept>  // std::runtime_error

CRTBP2D::CRTBP2D(double mu, Model::Formalism formalism, Model::IndicatorType indicator)
{
    setName("Planar CRTBP");

    t_         = 0.0;  /// Elapsed dimensionless CRTBP time
    param_.mu  = mu;
    formalism_ = formalism;
    indicator_ = indicator;

    switch (indicator_) {
        case Model::IndicatorType::NONE:
            // Equations of motion only.
            setNVar(4);
            setFunction(&Model::fun);
            break;

        case Model::IndicatorType::FLI:
        case Model::IndicatorType::LCI:
            // Equations of motion and variational equations.
            setNVar(8);
            setFunction(&Model::varFun);
            break;

        case Model::IndicatorType::RLI:
            throw std::runtime_error("RLI indicator is not yet implemented.");

        default:
            throw std::runtime_error("Unknown indicator type.");
    }
}

void CRTBP2D::inertialToCRTBP(const astro::State &state, double a2, double n)
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

astro::State CRTBP2D::crtbpToInertial(const double *y, double a2, double n) const noexcept
{
    const double mu = param_.mu;

    const double x  = y[0];
    const double yr = y[1];
    const double vx = y[2];
    const double vy = y[3];

    // Position relative to P1 in the rotating frame.
    const double xP1 = x + mu;
    const double yP1 = yr;

    // Inertial velocity relative to P1, expressed in rotating coordinates.
    const double vxInert = vx - yP1;
    const double vyInert = vy + xP1;

    // Rotation angle between the rotating and inertial frames.
    const double c = std::cos(t_);
    const double s = std::sin(t_);

    astro::State state{};

    // P1-centered inertial position [AU].
    state.r.x = a2 * (c * xP1 - s * yP1);
    state.r.y = a2 * (s * xP1 + c * yP1);
    state.r.z = 0.0;

    // P1-centered inertial velocity [AU/day].
    state.v.x = n * a2 * (c * vxInert - s * vyInert);
    state.v.y = n * a2 * (s * vxInert + c * vyInert);
    state.v.z = 0.0;

    return state;
}

void CRTBP2D::getInitialCondition(double a, double e)
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

void CRTBP2D::velocityToHamiltonian() noexcept
{
    double *y = getY();

    const double px = y[2] - y[1];
    const double py = y[3] + y[0];

    y[2] = px;
    y[3] = py;
}

void CRTBP2D::hamiltonianToNewtonian(double *y_out) const noexcept
{
    const double *y = getY();

    y_out[0] = y[0];
    y_out[1] = y[1];
    y_out[2] = y[2] + y[1];
    y_out[3] = y[3] - y[0];
}

void CRTBP2D::fun(double t, const double *y, double *dydt, void *par) const
{
    switch (formalism_) {
        case Formalism::NEWTONIAN:
            funNewtonian(t, y, dydt, par);
            break;

        case Formalism::HAMILTONIAN:
            funHamiltonian(t, y, dydt, par);
            break;

        default:
            throw std::runtime_error("Unknown CRTBP formalism.");
    }
}

void CRTBP2D::varFun(double t, const double *y, double *dydt, void *par) const
{
    switch (formalism_) {
        case Formalism::NEWTONIAN:
            varFunNewtonian(t, y, dydt, par);
            break;

        case Formalism::HAMILTONIAN:
            varFunHamiltonian(t, y, dydt, par);
            break;

        default:
            throw std::runtime_error("Unknown CRTBP formalism.");
    }
}

void CRTBP2D::funNewtonian(double t, const double *y, double *dydt, void *par) const
{
    (void)t;

    const double mu = param_.mu;

    const double r1 = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
    const double r2 = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));

    const double r1_3 = 1.0 / astro::cube(r1);
    const double r2_3 = 1.0 / astro::cube(r2);

    dydt[0] = y[2];
    dydt[1] = y[3];
    dydt[2] = 2.0 * y[3] + y[0] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;
    dydt[3] = -2.0 * y[2] + y[1] * (1.0 - (1.0 - mu) * r1_3 - mu * r2_3);
}

void CRTBP2D::funHamiltonian(double t, const double *y, double *dydt, void *par) const
{
    (void)t;
    (void)par;

    const double mu = param_.mu;

    const double r1 = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
    const double r2 = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));

    const double r1_3 = 1.0 / astro::cube(r1);
    const double r2_3 = 1.0 / astro::cube(r2);

    // Hamiltonian equations of motion.
    dydt[0] = y[2] + y[1];
    dydt[1] = y[3] - y[0];
    dydt[2] = y[3] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;
    dydt[3] = -y[2] - (1.0 - mu) * y[1] * r1_3 - mu * y[1] * r2_3;
}

void CRTBP2D::varFunNewtonian(double t, const double *y, double *dydt, void *par) const
{
    (void)t;

    const auto  *p  = static_cast<const Params *>(par);
    const double mu = p->mu;

    const double r1   = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
    const double r2   = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));
    const double r1_3 = 1.0 / astro::cube(r1);
    const double r2_3 = 1.0 / astro::cube(r2);
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

void CRTBP2D::varFunHamiltonian(double t, const double *y, double *dydt, void *par) const
{
    (void)t;
    (void)par;

    const double mu = param_.mu;

    const double r1   = std::sqrt(astro::sqr(y[0] + mu) + astro::sqr(y[1]));
    const double r2   = std::sqrt(astro::sqr(y[0] - 1.0 + mu) + astro::sqr(y[1]));
    const double r1_3 = 1.0 / astro::cube(r1);
    const double r2_3 = 1.0 / astro::cube(r2);
    const double r1_5 = r1_3 / astro::sqr(r1);
    const double r2_5 = r2_3 / astro::sqr(r2);

    // Hamiltonian equations of motion.
    dydt[0] = y[2] + y[1];
    dydt[1] = y[3] - y[0];
    dydt[2] = y[3] - (1.0 - mu) * (y[0] + mu) * r1_3 - mu * (y[0] - 1.0 + mu) * r2_3;
    dydt[3] = -y[2] - (1.0 - mu) * y[1] * r1_3 - mu * y[1] * r2_3;

    // Second derivatives of the Hamiltonian.
    const double H_xx = (1.0 - mu) * r1_3 - 3.0 * (1.0 - mu) * astro::sqr(y[0] + mu) * r1_5 + mu * r2_3 -
                        3.0 * mu * astro::sqr(y[0] - 1.0 + mu) * r2_5;
    const double H_xy = -3.0 * (1.0 - mu) * (y[0] + mu) * y[1] * r1_5 - 3.0 * mu * (y[0] - 1.0 + mu) * y[1] * r2_5;
    const double H_yy =
        (1.0 - mu) * r1_3 - 3.0 * (1.0 - mu) * astro::sqr(y[1]) * r1_5 + mu * r2_3 - 3.0 * mu * astro::sqr(y[1]) * r2_5;

    // Hamiltonian variational equations.
    dydt[4] = y[6] + y[5];
    dydt[5] = y[7] - y[4];
    dydt[6] = -H_xx * y[4] - H_xy * y[5] + y[7];
    dydt[7] = -H_xy * y[4] - H_yy * y[5] - y[6];
}

void CRTBP2D::printState(std::ostream &os, double t, const double *y) const
{
    os << std::setw(io::DATA_FIELD_WIDTH) << t;
    for (std::size_t i = 0; i < getNVar(); ++i) {
        os << std::setw(io::DATA_FIELD_WIDTH) << y[i];
    }

    os << '\n';
}