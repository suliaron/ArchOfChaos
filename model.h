#pragma once

#include <cstddef>  // std::size_t
#include <memory>   // std::make_unique, std::unique_ptr
#include <ostream>  // std::ostream
#include <string>   // std::string

/**
 * @brief Abstract base class for dynamical models.
 *
 * Defines the common interface for evaluating the equations of motion
 * and the corresponding variational equations of a dynamical system.
 */
class Model {
   public:
    /**
     * @brief Type of a model right-hand-side member function.
     */
    using rhs_t = void (Model::*)(double t, const double *y, double *dydt, void *par) const;

    /**
     * @brief Virtual destructor.
     */
    virtual ~Model() = default;

    /**
     * @brief Evaluates the currently selected right-hand side.
     *
     * Calls the right-hand-side function previously selected by
     * setFunction().
     *
     * @param t Current time.
     * @param y State vector.
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to model-specific parameters.
     */
    void f(double t, const double *y, double *dydt, void *par) const
    {
        (this->*f_)(t, y, dydt, par);
    }

    /**
     * @brief Selects the right-hand-side function.
     *
     * @param f Pointer to the right-hand-side member function.
     */
    void setFunction(rhs_t f) noexcept
    {
        f_ = f;
    }

    /**
     * @brief Sets the number of dynamical variables and allocates the state vector.
     *
     * If the number of variables changes, the state vector is reallocated
     * to match the new size. Its previous contents are lost.
     *
     * @param n_var Number of dynamical variables.
     */
    void setNVar(std::size_t n_var)
    {
        if (n_var != n_var_) {
            n_var_ = n_var;
            y_ = std::make_unique<double[]>(n_var_);
        }
    }
    /**
     * @brief Returns the number of dynamical variables.
     *
     * @return Number of variables.
     */
    std::size_t getNVar() const noexcept
    {
        return n_var_;
    }

    /**
     * @brief Returns the state vector.
     *
     * @return Pointer to the state vector.
     */
    double* getY() noexcept
    {
        return y_.get();
    }

    /**
     * @brief Returns the state vector for read-only access.
     *
     * @return Const pointer to the state vector.
     */
    const double* getY() const noexcept
    {
        return y_.get();
    }

    /**
     * @brief Sets the model name.
     *
     * @param name Model name.
     */
    void setName(const std::string &name)
    {
        name_ = name;
    }

    /**
     * @brief Returns the model name.
     *
     * @return Model name.
     */
    const std::string &getName() const noexcept
    {
        return name_;
    }

    /**
     * @brief Evaluates the equations of motion.
     *
     * @param t Current time.
     * @param y State vector.
     * @param dydt Time derivative of the state vector.
     * @param par Pointer to model-specific parameters.
     */
    virtual void fun(double t, const double *y, double *dydt, void *par) const = 0;

    /**
     * @brief Evaluates the equations of motion and variational equations.
     *
     * @param t Current time.
     * @param y State vector including the deviation vector.
     * @param dydt Time derivative of the state and deviation vectors.
     * @param par Pointer to model-specific parameters.
     */
    virtual void varfun(double t, const double *y, double *dydt, void *par) const = 0;

    /**
     * @brief Prints the current state of the model.
     *
     * @param os Output stream.
     * @param t Current time.
     * @param y State vector.
     */
    virtual void printState(std::ostream &os, double t, const double *y) const = 0;

   private:
    rhs_t       f_     = &Model::fun;
    std::size_t n_var_ = 0;
    std::string name_;
    std::unique_ptr<double[]> y_;
};