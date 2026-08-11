#pragma once

/**
 * @brief Base class for dynamical models.
 *
 * Defines the common interface for evaluating the equations of motion
 * and the corresponding variational equations of a dynamical system.
 */
class Model
{
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~Model() = default;

    /**
     * @brief Evaluates the equations of motion.
     *
     * @param t Current time.
     * @param y State vector.
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to model-specific parameters.
     */
    virtual void fun(double t, const double* y, double* dydt, void* par) const = 0;

    /**
     * @brief Evaluates the equations of motion and the variational equations.
     *
     * @param t Current time.
     * @param y State vector including the deviation vector.
     * @param dydt Time derivative of the state and deviation vectors.
     * @param par Pointer to model-specific parameters.
     */
    virtual void varfun(double t, const double* y, double* dydt, void* par) const = 0;
};