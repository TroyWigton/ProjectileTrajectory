/**
 * @file simulation.hpp
 * @brief Defines the Simulation class for running projectile motion scenarios.
 */

#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "types.hpp"
#include "integrators.hpp"
#include "derivative_functions.hpp"

/**
 * @struct ScenarioResult
 * @brief Container for the results of a single simulation run.
 */
struct ScenarioResult {
    double angle;       ///< Launch angle in degrees
    double velocity;    ///< Initial velocity in m/s
    double height;      ///< Initial height in m
    double deltaT;      ///< Time step used in s
    double distance;    ///< Total horizontal distance traveled in m
    long long step_count; ///< Number of integration steps performed
};

/**
 * @class Simulation
 * @brief Manages configuration and execution of projectile trajectories.
 * 
 * Holds the physics parameters and integrator configuration.
 * Can be run repeatedly with different launch angles.
 */
class Simulation {
public:
    /**
     * @brief Constructs a new Simulation configuration.
     * 
     * @param k_over_m The drag coefficient divided by mass.
     * @param integrator The numerical integration strategy to use.
     * @param v0 Initial velocity magnitude (m/s).
     * @param derivative The physics model (derivative function) to use.
     * @param g Gravitational acceleration (m/s^2).
     * @param h0 Initial height (m).
     */
    Simulation(double k_over_m,
               SystemIntegrator<State4D> integrator,
               double v0 = 100.0,
               DerivativeFuncPtr derivative = drag_deriv_v_squared,
               double g = 9.81,
               double h0 = 0.0);

    /**
     * @brief Runs a simulation for a specific launch angle.
     * 
     * @param angle_deg Launch angle in degrees (relative to horizontal).
     * @param dt Time step for the integrator.
     * @return ScenarioResult struct containing the simulation outcome.
     */
    ScenarioResult run(double angle_deg, double dt) const;

    /** @return Configured gravitational constant. */
    double gravity() const { return g_; }
    /** @return Configured drag/mass ratio. */
    double drag_ratio() const { return k_over_m_; }
    /** @return Configured initial velocity. */
    double initial_velocity() const { return v0_; }
    /** @return Configured initial height. */
    double initial_height() const { return h0_; }

private:
    const double g_;         ///< Gravitational acceleration
    const double k_over_m_;  ///< Drag coefficient / mass
    const double v0_;        ///< Initial velocity magnitude
    const double h0_;        ///< Initial height
    
    // Abstracted stepper function that knows how to advance the system
    // (Integrator + Derivative + Parameters all bound)
    const std::function<State4D(const State4D&, double, double)> stepper_;
};

#endif // SIMULATION_HPP
