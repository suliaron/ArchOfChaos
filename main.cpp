#include "command_line.h"    /**< Command-line option parsing and storage. */
#include "crtbp.h"           /**< CRTBP2D and CRTBP3D models. */
#include "chaos_indicator.h" /**< Chaos-indicator calculations (FLI, LCI). */
#include "grid.h"            /**< Grid iteration and orbital-element grid handling. */
#include "init_data.h"       /**< Initialization data and run configuration. */
#include "integrator.h"      /**< Adaptive numerical integration and step-size control. */
#include "io.h"              /**< Input/output helpers and reproducibility metadata. */
#include "math_utils.h"      /**< Mathematical and astronomical utility functions. */
#include "model.h"           /**< Model base class and model-related types. */
#include "orbit.h"           /**< Orbital-element/state transformations. */
#include "time_utils.h"      /**< astro::LogOutputSchedule. */

#include <algorithm>  /**< std::any_of, std::copy_n, std::max, std::min. */
#include <chrono>     /**< Provides time measurement utilities. */
#include <cmath>      /**< Mathematical functions (std::sqrt, std::abs). */
#include <cstddef>    /**< std::size_t. */
#include <exception>  /**< std::exception base class. */
#include <filesystem> /**< std::filesystem::path for portable path handling. */
#include <fstream>    /**< File stream classes (std::ifstream, std::ofstream). */
#include <iomanip>    /**< Output manipulators (std::fixed, std::setprecision, std::setw). */
#include <iostream>   /**< Standard input/output streams (std::cout, std::cerr). */
#include <limits>     /**< std::numeric_limits  */
#include <stdexcept>  /**< Standard exception classes (std::invalid_argument, std::runtime_error). */
#include <vector>     /**< std::vector container. */

namespace run {
    namespace fs = std::filesystem;

    constexpr double TIME_EPS = 1.0e-15;

    /**
     * @brief Integrates a single CRTBP orbit and writes its state evolution.
     *
     * Constructs the initial heliocentric inertial state, transforms it to the
     * dimensionless rotating CRTBP frame, and integrates the equations of motion.
     * Numerical integration tolerances are obtained from the initialization data.
     *
     * @param model CRTBP model.
     * @param init Initialization data.
     * @param out Output stream.
     * @param step Adaptive integration step-control parameters.
     * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
     *             [AU^3/day^2].
     * @param n Mean motion of the P1-P2 system [rad/day].
     * @param outputDtDimless Dimensionless output time interval.
     *
     * @throws std::runtime_error If a non-finite state is encountered.
     */
    void orbit(Model &model, const InitData &init, std::ostream &out, StepControl &step, double mu_13, double n,
               double outputDtDimless)
    {
        // Copy the input orbital elements.
        astro::OrbitalElements elements = init.getElements();
        // Calculate the physical integration duration for the current orbit.
        const double duration = init.calcIntegrationDuration(mu_13, elements.a);
        // Convert the physical integration duration to dimensionless CRTBP time.
        const double tDimless = CRTBP2D::toDimlessTime(duration, n);
        // Calculate the pericenter passage time from the selected
        // orbital-phase input.
        elements.tau = init.calc_tau(mu_13, elements.a);

        // Orbital elements -> heliocentric inertial Cartesian state.
        const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

        // Heliocentric inertial state -> dimensionless rotating CRTBP state.
        model.inertialToCRTBP(state, init.getA2(), n);

        // Initial conditions from EB Andromeda p. 5. (these are rotated by 180 degrees)
        model.getY()[0] = -1.06201;
        model.getY()[1] = 0.0;
        model.getY()[2] = 0.0;
        model.getY()[3] = 0.10851;

        // Calculate the Jacobi constant for the initial Newtonian state.
        model.initializeJacobiConstant();

#if 0  // Test CRTBP -> inertial transformation
    std::ofstream inertialOut("output_inertial.txt");
    if (!inertialOut) {
        throw std::runtime_error("Cannot open output_inertial.txt.");
    }
    inertialOut << "# t [day]            x [AU]             y [AU]"
                   "             vx [AU/day]         vy [AU/day]\n";
    inertialOut << std::scientific << std::setprecision(10) << std::showpos;
    // Transform the initial CRTBP state back to the
    // P1-centered inertial reference frame.
    const astro::State initialInertial = model.crtbpToInertial(model.getY(), init.getA2(), n);
    inertialOut << init.getT0() << ' ' << initialInertial.r.x << ' ' << initialInertial.r.y << ' '
                << initialInertial.v.x << ' ' << initialInertial.v.y << '\n';
#endif

        // Save the initial state using physical time [day].
        model.printState(out, init.getT0(), model.getY());

        // Convert the orbital state to Hamiltonian canonical variables if required.
        if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
            model.velocityToHamiltonian();
        }

        // Set the first output time (in dimensionless CRTBP time) according to the integration direction.
        const double directionSign     = integrationDirectionSign(init.getIntegrationDirection());
        double       nextOutputDimless = directionSign * outputDtDimless;

        while (directionSign * (tDimless - model.getT()) > TIME_EPS) {
            if (step.n_tst % 10 == 0) {
                if (!checkFinite(model.getY(), model.getNVar())) {
                    throw std::runtime_error("Non-finite state encountered during orbit integration.");
                }
            }
            const double targetTime = (init.getIntegrationDirection() == IntegrationDirection::FORWARD)
                                          ? std::min(nextOutputDimless, tDimless)
                                          : std::max(nextOutputDimless, tDimless);
            // Force the integrator to stop exactly at the next output
            // or final dimensionless time.
            limitStep(model.getT(), targetTime, step);

            ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

            ++step.n_int;
            ++step.n_tst;

            const bool outputTimeReached = directionSign * (model.getT() - nextOutputDimless) >= -TIME_EPS;
            const bool finalTimeReached  = directionSign * (model.getT() - tDimless) >= -TIME_EPS;
            if (outputTimeReached || finalTimeReached) {
                const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);

                if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
                    std::vector<double> y_out(model.getNVar());
                    // Hamiltonian -> Newtonian state.
                    model.hamiltonianToNewtonian(y_out.data());
                    model.calcJacobiConstant();
                    // Save the Newtonian CRTBP state.
                    model.printState(out, physicalTime, y_out.data());

#if 0  // Test CRTBP -> inertial transformation
                const astro::State inertialState = model.crtbpToInertial(y_out.data(), init.getA2(), n);
                inertialOut << physicalTime << ' ' << inertialState.r.x << ' ' << inertialState.r.y << ' '
                            << inertialState.v.x << ' ' << inertialState.v.y << '\n';
#endif
                }
                else {
                    model.calcJacobiConstant();
                    // Save the Newtonian CRTBP state.
                    model.printState(out, physicalTime, model.getY());

#if 0  // Test CRTBP -> inertial transformation
                const astro::State inertialState = model.crtbpToInertial(model.getY(), init.getA2(), n);
                inertialOut << physicalTime << ' ' << inertialState.r.x << ' ' << inertialState.r.y << ' '
                            << inertialState.v.x << ' ' << inertialState.v.y << '\n';
#endif
                }
                if (outputTimeReached) {
                    nextOutputDimless += directionSign * outputDtDimless;
                }
            }
        } /* while */
    } /* orbit */

    /**
     * @brief Integrates a single CRTBP orbit and computes a chaos indicator.
     *
     * Initializes the orbital state and deviation vector, integrates the
     * equations of motion and variational equations in dimensionless CRTBP time,
     * and writes the selected chaos indicator at logarithmically distributed
     * physical output times.
     *
     * The first output time and the number of output points per decade are
     * specified by InitData. Output times represent elapsed physical times
     * measured from the initial epoch.
     *
     * The final integration point is always written, even if it does not
     * coincide with a scheduled output time.
     *
     * Physical input times are specified in days, while numerical integration
     * is performed using dimensionless CRTBP time.
     *
     * @param model CRTBP model.
     * @param init Initialization data.
     * @param out Output stream.
     * @param step Adaptive integration step-control parameters.
     * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit [AU^3/day^2].
     * @param n Mean motion of the P1-P2 system [rad/day].
     *
     * @throws std::runtime_error If a non-finite state is encountered.
     */
    void indicator(Model &model, const InitData &init, std::ostream &out, StepControl &step, double mu_13, double n)
    {
        // Copy the input orbital elements.
        astro::OrbitalElements elements = init.getElements();
        // Calculate the physical integration duration for the current orbit.
        const double duration = init.calcIntegrationDuration(mu_13, elements.a);
        // Convert the physical integration duration to dimensionless CRTBP time.
        const double tDimless = CRTBP2D::toDimlessTime(duration, n);
        // Calculate the pericenter passage time from the selected
        // orbital-phase input.
        elements.tau = init.calc_tau(mu_13, elements.a);
        // Orbital elements -> heliocentric inertial Cartesian state.
        const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

        // Heliocentric inertial state -> dimensionless rotating CRTBP state.
        model.inertialToCRTBP(state, init.getA2(), n);

        // Convert the orbital state to Hamiltonian canonical variables if required.
        if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
            model.velocityToHamiltonian();
        }

        // Number of components in the deviation vector.
        const std::size_t nDeviation = model.getNVar() / 2;

        // Set the initial deviation vector.
        double *deviation = model.getY() + nDeviation;
        std::copy_n(init.getDy(), nDeviation, deviation);

        // Initial norm of the deviation vector.
        const double norm_0 = astro::norm(deviation, nDeviation);

        // ---------------------------------------------------------------------
        // Initialize the selected chaos indicator.
        // ---------------------------------------------------------------------
        double indicator_value = 0.0;
        switch (model.getIndicator()) {
            case Model::IndicatorType::FLI:
                // Initialize the running FLI value at the initial time.
                //
                // The initial value itself is not written because t = 0 cannot
                // be represented on a logarithmic time axis.
                indicator_value = chaos_indicator::computeFLI(deviation, nDeviation, indicator_value, norm_0);
                break;

            case Model::IndicatorType::LCI:
                // LCI is not defined at the initial time.
                break;

            case Model::IndicatorType::RLI:
                throw std::runtime_error("RLI indicator is not yet implemented.");

            case Model::IndicatorType::NONE:
                throw std::runtime_error("INDICATOR mode requires a chaos indicator.");

            default:
                throw std::runtime_error("Unknown chaos indicator.");
        }

        // ---------------------------------------------------------------------
        // Construct the logarithmic output schedule.
        // ---------------------------------------------------------------------
        astro::LogOutputSchedule outputSchedule(init.getOutputFirst(), init.getOutputPointsPerDecade());
        // Integration direction: +1 for forward and -1 for backward integration.
        const double directionSign = integrationDirectionSign(init.getIntegrationDirection());
        // First scheduled output time in dimensionless CRTBP units.
        double nextOutputDimless = directionSign * CRTBP2D::toDimlessTime(outputSchedule.getNextTime(), n);

        // ---------------------------------------------------------------------
        // Integrate the orbit and variational equations.
        // ---------------------------------------------------------------------
        while (directionSign * (tDimless - model.getT()) > TIME_EPS) {
            // Check the numerical state periodically.
            if (step.n_tst % 10 == 0) {
                if (!checkFinite(model.getY(), model.getNVar())) {
                    throw std::runtime_error("Non-finite state encountered during indicator integration.");
                }
            }
            // Stop exactly at the next scheduled output time or at the final
            // integration time, whichever comes first in the integration direction.
            const double targetTime = (init.getIntegrationDirection() == IntegrationDirection::FORWARD)
                                          ? std::min(nextOutputDimless, tDimless)
                                          : std::max(nextOutputDimless, tDimless);

            limitStep(model.getT(), targetTime, step);

            // Integrate the orbit and variational equations.
            ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

            ++step.n_int;
            ++step.n_tst;

            // Physical epoch corresponding to the current dimensionless
            // CRTBP integration time.
            const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);

            // -----------------------------------------------------------------
            // Update the chaos indicator.
            // -----------------------------------------------------------------

            switch (model.getIndicator()) {
                case Model::IndicatorType::FLI:
                    // FLI is a running maximum and must be updated after every
                    // accepted integration step.
                    indicator_value = chaos_indicator::computeFLI(deviation, nDeviation, indicator_value, norm_0);
                    break;

                case Model::IndicatorType::LCI:
                    // Calculate the LCI using physical time so that its unit
                    // remains 1/day.
                    indicator_value =
                        chaos_indicator::computeLCI(physicalTime, init.getT0(), deviation, init.getDy(), nDeviation);
                    break;

                default:
                    break;
            }

            // -----------------------------------------------------------------
            // Check whether an output point has been reached.
            // -----------------------------------------------------------------

            const bool outputTimeReached = directionSign * (model.getT() - nextOutputDimless) >= -TIME_EPS;
            const bool finalTimeReached  = directionSign * (model.getT() - tDimless) >= -TIME_EPS;
            if (outputTimeReached || finalTimeReached) {
                // Write the current physical epoch and chaos-indicator value.
                //
                // If the final time is also a scheduled output time, the value
                // is written only once because both conditions are handled here.
                out << std::setw(io::DATA_FIELD_WIDTH) << physicalTime << std::setw(io::DATA_FIELD_WIDTH)
                    << indicator_value << '\n';

                // Advance the logarithmic schedule only when a scheduled
                // output point has actually been reached.
                if (outputTimeReached) {
                    outputSchedule.advance();
                    nextOutputDimless = directionSign * CRTBP2D::toDimlessTime(outputSchedule.getNextTime(), n);
                }
            }
        }
    }

    /**
     * @brief Computes a chaos indicator over a multidimensional
     *        orbital-element grid.
     *
     * Integrates the CRTBP equations and variational equations for every
     * point of an arbitrary orbital-element grid and writes the final value
     * of the selected chaos indicator.
     *
     * Fixed orbital elements are obtained from InitData, while grid-controlled
     * elements are replaced by the current values supplied by GridIterator.
     *
     * The orbital phase may be specified either by tau or M, as a fixed input
     * value or as a grid axis. If M is a grid axis, the corresponding time of
     * pericenter passage is calculated from the current semimajor axis.
     *
     * If the integration duration is specified by nPeriods, the physical
     * duration is recalculated independently at every grid point from the
     * current semimajor axis of P3.
     *
     * Physical times are expressed in days, while numerical integration is
     * performed in dimensionless CRTBP time.
     *
     * @param model CRTBP model.
     * @param init Initialization data.
     * @param out Output stream.
     * @param mu_13 Gravitational parameter of the P1-P3 heliocentric orbit
     *              [AU^3/day^2].
     * @param n Mean motion of the P1-P2 system [rad/day].
     *
     * @throws std::runtime_error If the selected chaos indicator is invalid.
     */
    void grid(Model &model, const InitData &init, std::ostream &out, double mu_13, double n)
    {
        // Construct the multidimensional orbital-element grid iterator.
        GridIterator grid(init.getGridAxes());
        const double directionSign = integrationDirectionSign(init.getIntegrationDirection());
        // Number of components in the deviation vector.
        const std::size_t nDeviation = model.getNVar() / 2;

        // ---------------------------------------------------------------------
        // Check the selected chaos indicator.
        // ---------------------------------------------------------------------
        switch (model.getIndicator()) {
            case Model::IndicatorType::FLI:
            case Model::IndicatorType::LCI:
                break;

            case Model::IndicatorType::RLI:
                throw std::runtime_error("RLI indicator is not yet implemented.");

            case Model::IndicatorType::NONE:
                throw std::runtime_error("GRID mode requires a chaos indicator.");

            default:
                throw std::runtime_error("Unknown chaos indicator.");
        }

        // ---------------------------------------------------------------------
        // Iterate over all grid points.
        // ---------------------------------------------------------------------
        do {
            // Reset the dimensionless CRTBP integration time.
            model.setT(0.0);
            // -----------------------------------------------------------------
            // Construct the orbital elements for the current grid point.
            // -----------------------------------------------------------------
            // Start from the fixed orbital elements specified in the input file.
            astro::OrbitalElements elements = init.getElements();
            // Apply all grid-controlled orbital elements that are stored
            // directly in astro::OrbitalElements.
            grid.apply(elements);
            // Initialize the adaptive step control for the current grid point.
            // The initial step size depends on the current semimajor axis of P3.
            StepControl step = createStepControl(elements.a, mu_13, n, init.getIntegrationDirection());

            // -----------------------------------------------------------------
            // Resolve the orbital phase.
            // -----------------------------------------------------------------
            if (grid.hasAxis(OrbitalElement::MEAN_ANOMALY)) {
                // Mean anomaly is specified in degrees in the grid.
                const double meanAnomaly = astro::toRad(grid.getValue(OrbitalElement::MEAN_ANOMALY));

                // Keplerian mean motion of P3 [rad/day].
                const double n_3 = std::sqrt(mu_13 / astro::cube(elements.a));

                // M(t0) = n_3 * (t0 - tau)
                // therefore
                // tau = t0 - M(t0) / n_3.
                elements.tau = init.getT0() - meanAnomaly / n_3;
            }
            else if (grid.hasAxis(OrbitalElement::PERICENTER_TIME)) {
                // tau has already been assigned by grid.apply().
                // No additional conversion is required.
            }
            else {
                // The orbital phase is fixed in the input file.
                //
                // calc_tau() returns the specified tau directly or calculates
                // it from the fixed mean anomaly using the current semimajor axis.
                elements.tau = init.calc_tau(mu_13, elements.a);
            }

            // -----------------------------------------------------------------
            // Calculate the integration duration for the current grid point.
            // -----------------------------------------------------------------
            const double duration = init.calcIntegrationDuration(mu_13, elements.a);
            const double tDimless = CRTBP2D::toDimlessTime(duration, n);

            // -----------------------------------------------------------------
            // Construct the initial CRTBP state.
            // -----------------------------------------------------------------

            // Orbital elements -> heliocentric inertial Cartesian state.
            const astro::State state = astro::calcState(mu_13, init.getT0(), elements);

            // Heliocentric inertial state ->
            // dimensionless rotating CRTBP state.
            model.inertialToCRTBP(state, init.getA2(), n);

            // Convert only the orbital state to Hamiltonian canonical
            // variables if required.
            if (model.getFormalism() == Model::Formalism::HAMILTONIAN) {
                model.velocityToHamiltonian();
            }

            // -----------------------------------------------------------------
            // Initialize the deviation vector.
            // -----------------------------------------------------------------
            double *deviation = model.getY() + nDeviation;
            // Set the initial deviation vector exactly as specified
            // in the input file.
            std::copy_n(init.getDy(), nDeviation, deviation);
            // Initial deviation-vector norm for the current grid point.
            const double norm_0 = astro::norm(deviation, nDeviation);

            // -----------------------------------------------------------------
            // Initialize the selected chaos indicator.
            // -----------------------------------------------------------------
            double indicator_value = 0.0;

            if (model.getIndicator() == Model::IndicatorType::FLI) {
                // At the initial time:
                //
                // log(||delta(t0)|| / ||delta(t0)||) = 0.
                //
                // Calling computeFLI() here also validates the initial
                // deviation-vector norm.
                indicator_value = chaos_indicator::computeFLI(deviation, nDeviation, 0.0, norm_0);
            }

            // -----------------------------------------------------------------
            // Integrate the current grid point.
            // -----------------------------------------------------------------

            bool valid = true;

            while (directionSign * (tDimless - model.getT()) > TIME_EPS) {
                // Check the numerical state periodically.
                if (step.n_tst % 10 == 0) {
                    if (!checkFinite(model.getY(), model.getNVar())) {
                        valid = false;
                        break;
                    }
                }

                // Force the last integration step to end exactly at
                // the final dimensionless integration time.
                limitStep(model.getT(), tDimless, step);

                // Integrate the orbit and variational equations.
                ode_integrator::rkf54(model, model.getParams(), step, init.getRelTol(), init.getAbsTol());

                ++step.n_int;
                ++step.n_tst;

                // FLI is a running maximum and must therefore be updated
                // after every accepted integration step.
                if (model.getIndicator() == Model::IndicatorType::FLI) {
                    indicator_value = chaos_indicator::computeFLI(deviation, nDeviation, indicator_value, norm_0);
                }
            } /* while */

            // -----------------------------------------------------------------
            // Final indicator value.
            // -----------------------------------------------------------------

            if (valid) {
                // LCI only needs to be evaluated at the final integration time.
                if (model.getIndicator() == Model::IndicatorType::LCI) {
                    // Convert the final dimensionless CRTBP time to physical time [day].
                    const double physicalTime = init.getT0() + CRTBP2D::toPhysicalTime(model.getT(), n);
                    indicator_value =
                        chaos_indicator::computeLCI(physicalTime, init.getT0(), deviation, init.getDy(), nDeviation);
                }
            }
            else {
                indicator_value = std::numeric_limits<double>::quiet_NaN();

                std::cerr << "\nNon-finite state encountered at grid point:\n";
                grid.printCurrentPoint(std::cerr);
                std::cerr << "Proceeding to the next grid point.\n";
            }

            // -----------------------------------------------------------------
            // Write the current grid coordinates and the indicator value.
            // -----------------------------------------------------------------
            for (std::size_t axisIndex = 0; axisIndex < grid.getAxisCount(); ++axisIndex) {
                out << std::setw(io::DATA_FIELD_WIDTH) << grid.getValue(axisIndex);
            }
            out << std::setw(io::DATA_FIELD_WIDTH) << indicator_value << '\n';

            // Display the grid progress.
            grid.printProgress(std::cerr);
        } while (grid.next());

        std::cerr << '\n';
    }

    void execute(const InitData &init, const CommandLineOptions &opt)
    {
        const double mu    = init.getM2() / (init.getM1() + init.getM2());
        const double mu_12 = astro::sqr(astro::k) * (init.getM1() + init.getM2());
        const double mu_13 = astro::sqr(astro::k) * init.getM1();
        const double n     = std::sqrt(mu_12 / astro::cube(init.getA2()));

        Model *model;
        switch (init.getProblemType()) {
            case ProblemType::CRTBP2D:
                model = new CRTBP2D(mu, init.getFormalism(), init.getIndicator());
                break;
            case ProblemType::CRTBP3D:
                model = new CRTBP3D(mu, init.getFormalism(), init.getIndicator());
                break;
            default:
                throw std::runtime_error("Unknown problem type.");
        }

        // Open either the explicitly requested output file or the automatically
        // generated output file.
        std::ofstream fout;
        std::ostream *out = io::openOutputStream(opt, init, fout);

        // Write the complete reproducibility header:
        //  - program and run metadata,
        //  - exact copy of the input file,
        //  - description of the numerical output structure.
        io::writeOutputHeader(*out, fs::path(opt.input_path), init);
        io::configureNumericalOutput(*out);
        // Write the actual numerical table header from the same schema that is
        // used in the reproducibility header.
        io::writeTableHeader(*out, init);
        out->flush();

        // ---------------------------------------------------------------------
        // Execute the selected computation mode.
        // ---------------------------------------------------------------------
        switch (init.getRunMode()) {
            case RunMode::ORBIT: {
                // Dimensionless output time interval for ORBIT mode.
                const double outputDtDimless = CRTBP2D::toDimlessTime(init.getOutputDt(), n);
                const double a               = init.getElements().a;
                StepControl  step            = createStepControl(a, mu_13, n, init.getIntegrationDirection());
                orbit(*model, init, *out, step, mu_13, n, outputDtDimless);
                break;
            }

            case RunMode::INDICATOR: {
                const double a    = init.getElements().a;
                StepControl  step = createStepControl(a, mu_13, n, init.getIntegrationDirection());
                indicator(*model, init, *out, step, mu_13, n);
                break;
            }

            case RunMode::GRID:
                grid(*model, init, *out, mu_13, n);
                break;

            default:
                delete model;
                throw std::runtime_error("Unknown run mode.");
        }
        delete model;
    }

}  // namespace run

int main(int argc, char *argv[])
{
    try {
        CommandLineOptions opt;
        parseCommandLine(argc, argv, opt);

        if (opt.show_version) {
            printVersion();
            return EXIT_SUCCESS;
        }
        if (opt.show_help) {
            printHelp();
            return EXIT_SUCCESS;
        }

        const InitData init(opt.input_path);

        if (opt.verbose) {
            init.print(std::cout);
        }

        const auto start_time = std::chrono::steady_clock::now();

        run::execute(init, opt);

        const auto                          end_time     = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed_time = end_time - start_time;
        std::cout << "\nTotal runtime : " << elapsed_time.count() << " s\n";

        return EXIT_SUCCESS;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << '\n';
    }
    catch (...) {
        std::cerr << "Unknown error.\n";
    }

    return EXIT_FAILURE;
}
