#pragma once

#include "model.h"  // Model base class

/**
 * @brief Planar circular restricted three-body problem model.
 *
 * Implements the equations of motion and the corresponding variational
 * equations of the planar circular restricted three-body problem (CRTBP)
 * in the rotating reference frame.
 */
class CRTBP2D : public Model
{
public:
    /**
     * @brief Parameters of the planar circular restricted three-body problem.
     */
    struct CRTBP2DParams
    {
        double mu;  /**< Gravitational mass parameter. */
    };

    /**
     * @brief Evaluates the equations of motion of the planar CRTBP.
     *
     * @param t Current time.
     * @param y State vector (x, y, vx, vy).
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to a CRTBP2DParams structure.
     */
    void fun(double t, const double* y, double* dydt, void* par) const override;

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
    void varfun(double t, const double* y, double* dydt, void* par) const override;
};