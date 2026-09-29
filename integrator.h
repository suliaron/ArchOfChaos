#pragma once

#include "init_data.h"
#include "model.h"

#include <cstddef>
#include <cstdint>

/**
 * @brief Adaptive step-size control parameters.
 *
 * Stores the current, accepted, proposed, and allowed integration
 * step sizes used by the adaptive Runge-Kutta integrator.
 */
struct StepControl {
    double h     = 0.0;  ///< Current integration step size.
    double h_nxt = 0.0;  ///< Proposed step size for the next step.
    double h_did = 0.0;  ///< Accepted step size of the current step.
    double h_min = 0.0;  ///< Minimum allowed step size.

    std::uint32_t n_int = 0;  ///< Number of integration steps taken.
    std::uint32_t n_tst = 0;  ///< Number of integration-step tests.
};

namespace ode_integrator {

    /**
     * @brief Performs one adaptive Runge-Kutta-Fehlberg 5(4) integration step.
     *
     * Advances the state stored in @p model using an embedded seven-stage
     * Runge-Kutta-Fehlberg method. The local truncation error is estimated
     * from the embedded pair and used to adapt the integration step size.
     *
     * Internal work arrays are stored in static vectors and are resized only
     * when the number of model variables changes.
     *
     * @param model Model containing the state vector and equations of motion.
     * @param par Model-specific parameter structure.
     * @param step Adaptive step-size control parameters.
     * @param relTol Relative error tolerance.
     * @param absTol Absolute error tolerance.
     */
    void rkf54(Model &model, void *par, StepControl &step, double relTol, double absTol);

}  // namespace ode_integrator

/**
 * @brief Creates the initial adaptive integration step-control parameters.
 *
 * The initial step size is chosen as a fixed fraction of the initial
 * Keplerian orbital period of P3 and converted to dimensionless CRTBP time.
 * Its sign is determined by the integration direction.
 *
 * @param a Semimajor axis of P3 [AU].
 * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
 *              [AU^3/day^2].
 * @param n Mean motion of the P1-P2 system [rad/day].
 * @param direction Direction of numerical time integration.
 *
 * @return Initialized step-control parameters.
 */
StepControl createStepControl(double a, double mu_13, double n, IntegrationDirection direction);

/**
 * @brief Limits the current integration step to a target time.
 *
 * Ensures that the current integration step does not pass the specified
 * target time. The function works for both forward and backward integration.
 *
 * @param t Current dimensionless integration time.
 * @param targetTime Target dimensionless integration time.
 * @param step Adaptive step-size control parameters.
 */
void limitStep(double t, double targetTime, StepControl &step);

/**
 * @brief Checks whether all elements of a state vector are finite.
 *
 * @param y State vector.
 * @param n Number of elements in the state vector.
 *
 * @return true if all elements are finite, false otherwise.
 */
bool checkFinite(const double *y, std::size_t n);
