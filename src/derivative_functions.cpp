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

void variable_drag_deriv(const State& s, double t, State& deriv, double g, double base_k_over_m) {
    // 1. Calculate current speed from State
    double v_sq = s[X_VEL]*s[X_VEL] + s[Y_VEL]*s[Y_VEL];
    double v = std::sqrt(v_sq);
    double mach = v / 343.0; // Approx speed of sound at sea level

    // 2. Define Physics Logic: Drag Multiplier vs Mach Number
    // Model: Rise from 1.0x at Mach 0.8 to 2.5x at Mach 1.0 (Transonic)
    //        Then decay back to 1.0x at Mach 2.0 (Supersonic)
    double drag_multiplier = 1.0;

    if (mach >= 0.8 && mach < 1.0) {
        // Linear rise: (mach - 0.8) / 0.2 gives 0..1 fraction
        drag_multiplier = 1.0 + ((mach - 0.8) / 0.2) * 1.5;
    } else if (mach >= 1.0 && mach < 2.0) {
        // Linear decay: (mach - 1.0) / 1.0 gives 0..1 fraction
        drag_multiplier = 2.5 - ((mach - 1.0) / 1.0) * 1.5;
    } else if (mach >= 2.0) {
        // Settled back to base at high mach (per request)
        // Note: Realistically it might stay higher, but we fade to base here
        drag_multiplier = 1.0;
    }
    
    double current_k_over_m = base_k_over_m * drag_multiplier;

    // 3. Apply force
    deriv[X_POS] = s[X_VEL];
    deriv[Y_POS] = s[Y_VEL];
    
    // Apply drag acceleration
    deriv[X_VEL] = -current_k_over_m * v * s[X_VEL];
    deriv[Y_VEL] = -g - current_k_over_m * v * s[Y_VEL];
}