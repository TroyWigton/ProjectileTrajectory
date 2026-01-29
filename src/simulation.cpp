/**
 * @file simulation.cpp
 * @brief Core logic for running individual projectile trajectory simulations.
 *
 * This file implements the Simulation class, which coordinates the setup and
 * execution of a single projectile flight. It binds the specific integration method
 * (stepper) with the physical derivative model and system constants.
 *
 * Key functionalities:
 * - Initialization of the simulation environment (gravity, drag, wind, etc.).
 * - Execution loop that steps the simulation until ground impact (y < 0).
 * - Linear interpolation to find the precise landing coordinate between time steps.
 */
#include "../include/simulation.hpp"
#include "../include/constants.hpp"
#include <cmath>

Simulation::Simulation(double k_over_m,
                       SystemIntegrator integrator,
                       double v0,
                       DerivativeFuncPtr derivative,
                       double g,
                       double h0)
    : g_(g),
      k_over_m_(k_over_m),
      v0_(v0),
      h0_(h0) {
    // initialize stepper by binding parameters and derivative to the integrator
    stepper_ = [integrator, derivative, g, k_over_m](const State& s, double t, double dt) -> State {
        // Lambda to match SystemDerivative signature
        auto bound_deriv = [derivative, g, k_over_m](const State& state, double time, State& out_d) {
             derivative(state, time, out_d, g, k_over_m); 
        };
        return integrator(s, t, dt, bound_deriv);
    };
}

ScenarioResult Simulation::run(double angle_deg, double dt) const {
    const double angle_rad = angle_deg * M_PI / 180.0;
    State state = {0.0, h0_, v0_ * std::cos(angle_rad), v0_ * std::sin(angle_rad)};
    double t = 0.0;
    double last_x = 0.0;
    double last_y = h0_;
    long long steps = 0;

    while (state[Y_POS] >= 0.0) {
        last_x = state[X_POS];
        last_y = state[Y_POS];
        state = stepper_(state, t, dt);
        t += dt;
        steps++;
    }

    const double denominator = state[Y_POS] - last_y;
    const double impact_x = std::fabs(denominator) > 1e-12
                                ? last_x + (-last_y) * (state[X_POS] - last_x) / denominator
                                : state[X_POS];

    return {angle_deg, v0_, h0_, dt, impact_x, steps};
}
