#pragma once

#include <cstddef>  // std::size_t
#include <limits>   // std::numeric_limits
#include <memory>   // std::make_unique, std::unique_ptr
#include <ostream>  // std::ostream
#include <string>   // std::string

// forward declaration:
namespace astro {
    struct State;
}

/**
 * @brief Abstract base class for dynamical models.
 *
 * Defines the common interface for evaluating the equations of motion
 * and the corresponding variational equations of a dynamical system.
 */
class Model {
   public:
    /**
     * @brief Specifies the mathematical formulation of the equations of motion.
     */
    enum class Formalism {
        NEWTONIAN,  /**< Position-velocity formulation. */
        HAMILTONIAN /**< Canonical Hamiltonian formulation. */
    };

    /**
     * @brief Specifies the chaos indicator to be computed.
     */
    enum class IndicatorType {
        NONE, /**< No chaos indicator is computed. */
        FLI,  /**< Fast Lyapunov Indicator. */
        LCI,  /**< Lyapunov Characteristic Indicator. */
        RLI   /**< Relative Lyapunov Indicator. */
    };

    /**
     * @brief Returns the name of a mathematical formulation.
     *
     * @param formalism Mathematical formulation.
     * @return Name of the formulation.
     */
    static const char *formalismToString(Formalism formalism) noexcept;

    /**
     * @brief Returns the name of a chaos indicator.
     *
     * @param indicator Chaos indicator type.
     * @return Name of the indicator.
     */
    static const char *indicatorTypeToString(IndicatorType indicator) noexcept;

    /**
     * @brief Type of a model right-hand-side member function.
     */
    using rhs_t = void (Model::*)(double t, const double *y, double *dydt, void *par) const;

    Model(Formalism formalism, IndicatorType indicator) :
        formalism_(formalism),
        indicator_(indicator)
    {
    }

    /**
     * @brief Virtual destructor.
     */
    virtual ~Model() = default;

    /**
     * @brief Returns the selected mathematical formulation.
     *
     * @return Current mathematical formulation.
     */
    Model::Formalism getFormalism() const noexcept
    {
        return formalism_;
    }

    /**
     * @brief Returns the selected chaos indicator.
     *
     * @return Current chaos indicator type.
     */
    Model::IndicatorType getIndicator() const noexcept
    {
        return indicator_;
    }

    /**
     * @brief Calculates the Jacobi constant of the system.
     *
     * @details This is a virtual fallback method. Derived classes should override
     * this function to provide actual computation. If a derived model does not
     * support or define the Jacobi constant, this default implementation returns NaN.
     *
     * @return The calculated Jacobi constant, or `std::numeric_limits<double>::quiet_NaN()`
     *         if not applicable/defined for the specific model.
     */
    virtual double calcJacobiConstant() const noexcept
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * @brief Initializes the Jacobi constant.
     *
     * Computes the Jacobi constant of the current state and stores it
     * both as the initial and current value.
     */
    void initializeJacobiConstant() noexcept
    {
        cj_ = calcJacobiConstant();
        c0_ = cj_;
    }

    /**
     * @brief Updates the current Jacobi constant.
     */
    void updateJacobiConstant() noexcept
    {
        cj_ = calcJacobiConstant();
    }

    /**
     * @brief Returns the model-specific parameter block.
     *
     * The returned pointer can be passed directly to the numerical
     * integrator and the model right-hand-side functions.
     *
     * @return Pointer to the model-specific parameters.
     */
    virtual void *getParams() noexcept = 0;

    /**
     * @brief Converts a heliocentric inertial Cartesian state to the model state.
     *
     * @param state Heliocentric inertial Cartesian state.
     * @param a2 Constant distance between the primary bodies [AU].
     * @param n Mean motion of the primary bodies [rad/day].
     */
    virtual void inertialToCRTBP(const astro::State &state, double a2, double n) = 0;

    /**
     * @brief Converts the model state from velocity variables to Hamiltonian
     *        canonical variables.
     */
    virtual void velocityToHamiltonian() noexcept = 0;

    /**
     * @brief Converts the current Hamiltonian state to Newtonian
     *        position-velocity variables.
     *
     * @param y_out Output Newtonian state vector.
     */
    virtual void hamiltonianToNewtonian(double *y_out) const noexcept = 0;

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
            y_     = std::make_unique<double[]>(n_var_);
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
     * @brief Sets the current model time.
     *
     * @param t Current time.
     */
    void setT(double t) noexcept
    {
        t_ = t;
    }

    /**
     * @brief Returns the current model time.
     *
     * @return Current time.
     */
    double getT() const noexcept
    {
        return t_;
    }

    /**
     * @brief Returns the state vector.
     *
     * @return Pointer to the state vector.
     */
    double *getY() noexcept
    {
        return y_.get();
    }

    /**
     * @brief Returns the state vector for read-only access.
     *
     * @return Const pointer to the state vector.
     */
    const double *getY() const noexcept
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
     * @brief Returns the initial Jacobi constant.
     *
     * @return Initial Jacobi constant.
     */
    double getC0() const noexcept
    {
        return c0_;
    }

    /**
     * @brief Returns the current Jacobi constant.
     *
     * @return Current Jacobi constant.
     */
    double getCJ() const noexcept
    {
        return cj_;
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
    virtual void varFun(double t, const double *y, double *dydt, void *par) const = 0;

    /**
     * @brief Prints the current state of the model.
     *
     * @param os Output stream.
     * @param t Current time.
     * @param y State vector.
     */
    virtual void printState(std::ostream &os, double t, const double *y) const = 0;

   protected:
    /**
     * @brief Sets the initial Jacobi constant.
     *
     * @param value Initial Jacobi constant.
     */
    void setC0(double value) noexcept
    {
        c0_ = value;
    }

    /**
     * @brief Sets the current Jacobi constant.
     *
     * @param value Current Jacobi constant.
     */
    void setCJ(double value) noexcept
    {
        cj_ = value;
    }

    Formalism                 formalism_ = Formalism::NEWTONIAN;
    IndicatorType             indicator_ = IndicatorType::NONE;
    double                    t_         = 0.0;
    rhs_t                     f_         = &Model::fun;
    std::size_t               n_var_     = 0;
    std::string               name_;
    std::unique_ptr<double[]> y_;

   private:
    double c0_ = std::numeric_limits<double>::quiet_NaN();
    double cj_ = std::numeric_limits<double>::quiet_NaN();
};