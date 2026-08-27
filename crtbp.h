#pragma once

#include "astro_types.h"
#include "model.h"  // Model base class

#include <ostream>  // std::ostream
/**
 * @brief Planar circular restricted three-body problem model.
 *
 * Implements the equations of motion and the corresponding variational
 * equations of the planar circular restricted three-body problem (CRTBP)
 * in the rotating reference frame.
 */
class CRTBP2D : public Model {
   public:
    /**
     * @brief Parameters of the planar circular restricted three-body problem.
     */
    struct Params {
        double mu; /**< CRTBP mass parameter, mu = m2 / (m1 + m2). */
    };

    /**
     * @brief Constructs a planar CRTBP model.
     *
     * Initializes the planar circular restricted three-body problem with the
     * specified mass parameter, mathematical formalism, and chaos indicator.
     *
     * The internal model time is initialized to zero and represents elapsed
     * dimensionless CRTBP time.
     *
     * @param mu CRTBP mass parameter.
     * @param formalism Mathematical formulation of the equations of motion.
     * @param indicator Chaos indicator to be computed.
     *
     * @throws std::runtime_error If the selected indicator is not implemented
     *         or is unknown.
     */
    CRTBP2D(double mu, Model::Formalism formalism, Model::IndicatorType indicator);

    /**
     * @brief Returns the model parameters.
     *
     * @return Pointer to the CRTBP parameters.
     */
    Params *getParams() noexcept
    {
        return &param_;
    }

    /**
     * @brief Converts a heliocentric inertial state to the normalized
     *        barycentric rotating CRTBP state.
     *
     * Converts a Cartesian state expressed in the heliocentric inertial
     * reference frame to the normalized barycentric rotating reference frame
     * used by the planar circular restricted three-body problem.
     *
     * The inertial position is assumed to be expressed in AU and the velocity
     * in AU/day. The CRTBP coordinates are dimensionless: distances are scaled
     * by the constant separation @p a2 of the two primary bodies, and time is
     * scaled by the inverse mean motion @p n.
     *
     * At the transformation epoch, the inertial and rotating coordinate axes
     * are assumed to be aligned, with the secondary body located on the
     * positive x-axis.
     *
     * The CRTBP mass parameter and the output state vector are taken directly
     * from the model.
     *
     * The resulting model state vector is
     *
     *     y = (x, y, vx, vy),
     *
     * where the position is barycentric and the velocity is measured in the
     * rotating reference frame.
     *
     * @param state Heliocentric inertial Cartesian state.
     * @param a2 Constant distance between the two primary bodies [AU].
     * @param n Mean motion of the primary bodies [rad/day].
     */
    void inertialToCRTBP(const astro::State &state, double a2, double n);

    /**
     * @brief Converts a Newtonian CRTBP state to the P1-centered inertial frame.
     *
     * Transforms a dimensionless barycentric rotating CRTBP state
     *
     *     y = (x, y, vx, vy)
     *
     * expressed in position-velocity variables to a heliocentric inertial
     * Cartesian state relative to the primary body P1.
     *
     * The rotating and inertial coordinate axes are assumed to be aligned at
     * dimensionless time t = 0. The current dimensionless model time is used
     * as the rotation angle between the two frames.
     *
     * @param y Newtonian CRTBP state vector (x, y, vx, vy).
     * @param a2 Constant distance between the two primary bodies [AU].
     * @param n Mean motion of the primary bodies [rad/day].
     *
     * @return P1-centered inertial Cartesian state with position in AU and
     *         velocity in AU/day.
     */
    astro::State crtbpToInertial(const double *y, double a2, double n) const noexcept;

    /**
     * @brief Sets the initial state for a special planar CRTBP configuration.
     *
     * Initializes the model state vector for a particle starting at the
     * periapsis of a planar elliptic orbit around the primary body. The
     * initial position lies on the positive x-axis, and the initial velocity
     * is perpendicular to the x-axis.
     *
     * The orbital elements are assumed to satisfy i = 0, omega = 0,
     * Omega = 0, and the initial time corresponds to the periapsis passage.
     * The semimajor axis @p a is expressed in normalized CRTBP distance units.
     *
     * The CRTBP mass parameter and the state vector are taken directly from
     * the model.
     *
     * @param a Semimajor axis of the particle orbit in normalized CRTBP units.
     * @param e Eccentricity of the particle orbit.
     */
    void getInitialCondition(double a, double e);

    /**
     * @brief Converts the model state to Hamiltonian canonical variables.
     *
     * Converts the normalized rotating CRTBP state from position-velocity
     * variables
     *
     *     y = (x, y, vx, vy)
     *
     * to Hamiltonian canonical variables
     *
     *     y = (x, y, px, py),
     *
     * using
     *
     *     px = vx - y,
     *     py = vy + x.
     *
     * The transformation is performed in place on the model state vector.
     * The position coordinates remain unchanged.
     */
    void velocityToHamiltonian() noexcept;

    /**
     * @brief Converts the Hamiltonian state to position-velocity variables.
     *
     * Converts the current normalized CRTBP state from canonical Hamiltonian
     * variables
     *
     *     y = (x, y, px, py)
     *
     * to position-velocity variables
     *
     *     y_out = (x, y, vx, vy),
     *
     * using
     *
     *     vx = px + y,
     *     vy = py - x.
     *
     * The internal model state is not modified. The converted state is written
     * to the output array @p y_out.
     *
     * @param y_out Output state vector (x, y, vx, vy).
     */
    void hamiltonianToNewtonian(double *y_out) const noexcept;

    /**
     * @brief Converts a physical time interval to dimensionless CRTBP time.
     *
     * The dimensionless CRTBP time interval is defined as
     *
     *     dt_dimless = n * dt_day,
     *
     * where @p n is the mean motion of the two primary bodies.
     *
     * @param dtDay Physical time interval [day].
     * @param n Mean motion of the primary bodies [rad/day].
     *
     * @return Time interval in dimensionless CRTBP units.
     */
    static double toDimlessTime(double dtDay, double n) noexcept
    {
        return n * dtDay;
    }

    /**
     * @brief Converts a dimensionless CRTBP time interval to physical time.
     *
     * The physical time interval is obtained from
     *
     *     dt_day = dt_dimless / n,
     *
     * where @p n is the mean motion of the two primary bodies.
     *
     * @param dtDimless Time interval in dimensionless CRTBP units.
     * @param n Mean motion of the primary bodies [rad/day].
     *
     * @return Physical time interval [day].
     */
    static double toPhysicalTime(double dtDimless, double n) noexcept
    {
        return dtDimless / n;
    }

    void printState(std::ostream &os, double t, const double *y) const override;

   private:
    /**
     * @brief Evaluates the equations of motion in the position-velocity formulation.
     *
     * The state vector is
     *
     *     y = (x, y, vx, vy).
     *
     * @param t Current time.
     * @param y State vector.
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to model parameters.
     */
    void funNewtonian(double t, const double *y, double *dydt, void *par) const;

    /**
     * @brief Evaluates the equations of motion in Hamiltonian canonical variables.
     *
     * The state vector is
     *
     *     y = (x, y, px, py).
     *
     * The normalized planar CRTBP equations are evaluated using
     * Omega = 1, GM1 = 1 - mu, and GM2 = mu.
     *
     * @param t Current time.
     * @param y State vector in canonical variables.
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to model parameters.
     */
    void funHamiltonian(double t, const double *y, double *dydt, void *par) const;

    /**
     * @brief Evaluates the equations of motion using the selected formalism.
     *
     * Dispatches the evaluation to either the Newtonian position-velocity
     * formulation or the Hamiltonian canonical formulation according to the
     * currently selected formalism.
     *
     * @param t Current time.
     * @param y State vector.
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to model-specific parameters.
     */
    void fun(double t, const double *y, double *dydt, void *par) const override;

    /**
     * @brief Evaluates the variational equations in the position-velocity formulation.
     *
     * The state vector contains the orbit and the corresponding deviation vector:
     *
     *     y = (x, y, vx, vy, dx, dy, dvx, dvy).
     *
     * @param t Current time.
     * @param y State vector and deviation vector.
     * @param dydt Time derivative of the state and deviation vectors.
     * @param par Pointer to model-specific parameters.
     */
    void varFunNewtonian(double t, const double *y, double *dydt, void *par) const;

    /**
     * @brief Evaluates the variational equations in Hamiltonian canonical variables.
     *
     * The state vector contains the orbit and the corresponding canonical
     * deviation vector:
     *
     *     y = (x, y, px, py, dx, dy, dpx, dpy).
     *
     * @param t Current time.
     * @param y State vector and deviation vector.
     * @param dydt Time derivative of the state and deviation vectors.
     * @param par Pointer to model-specific parameters.
     */
    void varFunHamiltonian(double t, const double *y, double *dydt, void *par) const;

    /**
     * @brief Evaluates the equations of motion and variational equations
     *        using the selected formalism.
     *
     * Dispatches the evaluation to either the Newtonian position-velocity
     * formulation or the Hamiltonian canonical formulation according to the
     * currently selected formalism.
     *
     * @param t Current time.
     * @param y State vector including the deviation vector.
     * @param dydt Time derivative of the state and deviation vectors.
     * @param par Pointer to model-specific parameters.
     */
    void varFun(double t, const double *y, double *dydt, void *par) const override;

    Params param_;
};