#include "integrator.h"

#include "crtbp.h"
#include "math_utils.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace ode_integrator {

    void rkf54(Model &model, void *par, StepControl &step, double relTol, double absTol)
    {
        static const double Bi[] = {17.0 / 192.0, 0.0, 64.0 / 231.0, 2187.0 / 8960.0, 2875.0 / 8448.0, 1.0 / 20.0, 0.0};
        static const double Ci[] = {0.0, 1.0 / 8.0, 1.0 / 4.0, 4.0 / 9.0, 4.0 / 5.0, 1.0, 1.0};
        static const double Aij[][6] = {
            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {1.0 / 8.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {196.0 / 729.0, -320.0 / 729.0, 448.0 / 729.0, 0.0, 0.0, 0.0},
            {836.0 / 2875.0, 64.0 / 575.0, -13376.0 / 20125.0, 21384.0 / 20125.0, 0.0, 0.0},
            {-73.0 / 48.0, 0.0, 1312.0 / 231.0, -2025.0 / 448.0, 2875.0 / 2112.0, 0.0},
            {17.0 / 192.0, 0.0, 64.0 / 231.0, 2187.0 / 8960.0, 2875.0 / 8448.0, 1.0 / 20.0}};

        const std::size_t n_var = model.getNVar();
        double           *y_in  = model.getY();

        static std::size_t         allocated_n_var = 0;
        static std::vector<double> dy;
        static std::vector<double> y;
        static std::vector<double> y_out;

        if (n_var != allocated_n_var) {
            dy.resize(7 * n_var);
            y.resize(n_var);
            y_out.resize(n_var);

            allocated_n_var = n_var;
        }

        double t0    = model.getT();
        double temax = 0.0;

        model.f(t0, y_in, dy.data(), par);
        do {
            temax = 0.0;

            for (std::size_t k = 1; k < 7; ++k) {
                const double t = t0 + Ci[k] * step.h;

                for (std::size_t n = 0; n < n_var; ++n) {
                    y[n] = y_in[n];

                    for (std::size_t l = 0; l < k; ++l) {
                        y[n] += step.h * Aij[k][l] * dy[l * n_var + n];
                    }
                }

                model.f(t, y.data(), dy.data() + k * n_var, par);
            }

            for (std::size_t n = 0; n < n_var; ++n) {
                y_out[n] = y_in[n];

                for (std::size_t k = 0; k < 7; ++k) {
                    y_out[n] += step.h * Bi[k] * dy[k * n_var + n];
                }

                const double err = std::abs(step.h) * std::fabs(dy[5 * n_var + n] - dy[6 * n_var + n]) / 60.0;
                const double tol = absTol + relTol * std::max(std::fabs(y_in[n]), std::fabs(y_out[n]));
                if (err / tol > temax) {
                    temax = err / tol;
                }
            }
            step.h_did = step.h;
            step.h = 0.9 * step.h_did * std::pow(1.0 / temax, 1.0 / 5.0);
        } while (temax > 1.0);
        step.h_nxt = step.h;

        // Advance to the end of the accepted integration step.
        model.setT(t0 + step.h_did);
        // Copy the accepted solution to the model state vector.
        std::copy_n(y_out.data(), n_var, model.getY());
    }
}  // namespace ode_integrator

StepControl createStepControl(double a, double mu_13, double n, IntegrationDirection direction)
{
    constexpr double INITIAL_STEPS_PER_PERIOD = 100.0;

    // Keplerian mean motion of P3 [rad/day].
    const double n_3 = std::sqrt(mu_13 / astro::cube(a));
    // Initial Keplerian orbital period of P3 [day].
    const double period = (2.0 * astro::pi) / n_3;
    // Convert the orbital period to dimensionless CRTBP time.
    const double periodDimless = CRTBP2D::toDimlessTime(period, n);
    // Initial adaptive integration step.
    const double h0 = periodDimless / INITIAL_STEPS_PER_PERIOD;

    StepControl step;
    step.h = integrationDirectionSign(direction) * h0;

    step.h_nxt = step.h;
    step.h_did = 0.0;
    step.n_tst = 0;
    step.n_int = 0;

    return step;
}

void limitStep(double t, double targetTime, StepControl &step)
{
    const double remaining = targetTime - t;

    if (std::abs(step.h) > std::abs(remaining)) {
        step.h = remaining;
    }
}

bool checkFinite(const double *y, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(y[i])) {
            return false;
        }
    }

    return true;
}
