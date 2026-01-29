/**
 * @file derivative_functions.cpp
 * @brief Implementation of physics models for projectile motion derivatives.
 *
 * This file contains the state derivative functions that describe the physical
 * equations of motion for the projectile system. These functions are used by
 * the numerical integrators to evolve the system state over time.
 *
 * Supported Models:
 * - v_squared_drag: Standard quadratic air resistance (F_drag ~ v^2).
 * - linear_drag: Linear air resistance (F_drag ~ v), useful for Stokes' law regime.
 * - no_drag: Vacuum trajectory, theoretical baseline.
 */
#include "../include/derivative_functions.hpp"
#include "../include/types.hpp"
#include <cmath>

void drag_deriv_v_squared(const State& s, double t, State& deriv, double g, double k_over_m) {
    const double v = std::sqrt(s[X_VEL]*s[X_VEL] + s[Y_VEL]*s[Y_VEL]);
    deriv[X_POS] = s[X_VEL];  // dx/dt = vx
    deriv[Y_POS] = s[Y_VEL];  // dy/dt = vy
    deriv[X_VEL] = -k_over_m * v * s[X_VEL];  // dvx/dt = -k/m * v * vx
    deriv[Y_VEL] = -g - k_over_m * v * s[Y_VEL];  // dvy/dt = -g -k/m * v * vy
}

void drag_deriv_linear(const State& s, double t, State& deriv, double g, double k_over_m) {
    deriv[X_POS] = s[X_VEL];  // dx/dt = vx
    deriv[Y_POS] = s[Y_VEL];  // dy/dt = vy
    deriv[X_VEL] = -k_over_m * s[X_VEL];  // dvx/dt = -k/m * vx
    deriv[Y_VEL] = -g - k_over_m * s[Y_VEL];  // dvy/dt = -g -k/m * vy
}

void no_drag_deriv(const State& s, double t, State& deriv, double g, double) {
    deriv[X_POS] = s[X_VEL];  // dx/dt = vx
    deriv[Y_POS] = s[Y_VEL];  // dy/dt = vy
    deriv[X_VEL] = 0.0;       // dvx/dt = 0
    deriv[Y_VEL] = -g;        // dvy/dt = -g
}