#pragma once

#include "astro_types.h"
#include "model.h"  // Model base class


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
    struct CRTBP2DParams {
        double mu; /**< Gravitational mass parameter. */
    };

    /**
     * @brief Constructs a planar CRTBP model.
     *
     * @param mu Gravitational mass parameter.
     */
    explicit CRTBP2D(double mu);

    /**
     * @brief Returns the model parameters.
     *
     * @return Pointer to the CRTBP parameters.
     */
    CRTBP2DParams *getParams() noexcept
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
    void InertialToCRTBP(const astro::State& state, double a2, double n);

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
    void GetInitialCondition(double a, double e);

    /**
     * @brief Evaluates the equations of motion of the planar CRTBP.
     *
     * @param t Current time.
     * @param y State vector (x, y, vx, vy).
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to a CRTBP2DParams structure.
     */
    void fun(double t, const double *y, double *dydt, void *par) const override;

    /**
     * @brief Evaluates the equations of motion and variational equations
     *        of the planar CRTBP.
     *
     * @param t Current time.
     * @param y State vector and deviation vector
     *          (x, y, vx, vy, dx, dy, dvx, dvy).
     * @param dydt Time derivative of the state and deviation vectors.
     * @param par Pointer to a CRTBP2DParams structure.
     */
    void varfun(double t, const double *y, double *dydt, void *par) const override;

    void printState(std::ostream &os, double t, const double *y) const override;

   private:
    CRTBP2DParams param_;
};