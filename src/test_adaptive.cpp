/**
 * @file test_adaptive.cpp
 * @brief Demonstration of adaptive Runge-Kutta 4/5 integration.
 *
 * This program runs a simulation of a golf ball trajectory using the RK45 adaptive
 * stepper. It demonstrates how to use the error estimate returned by the integrator
 * to dynamically adjust the time step (dt) during the simulation.
 *
 * This allows the simulation to take large steps when the physics are smooth,
 * and automatically slow down (take small steps) when complex phenomena occur
 * (like the golf ball drag crisis/transition).
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/integrators.hpp"

// Simple simulation runner for adaptive steps
void run_adaptive_simulation(double tolerance) {
    using namespace StateIndex4D;
    
    // Initial Conditions
    const double g = GRAVITY_EARTH;
    const double v0 = 70.0; // 70 m/s (~156 mph) - typical driver speed
    const double angle_deg = 15.0;
    const double angle_rad = angle_deg * M_PI / 180.0;
    const double h0 = 0.0;
    
    // Golf ball physics constants are handled inside the derivative function
    // But we need to pass a dummy dummy k/m or v_critical if the signature requires it.
    // golf_ball_drag_deriv signature: (s, t, deriv, g, v_critical)
    // We'll hardcode v_critical inside the loop or pass 0 if it's unused (the function uses shared constants)
    const double dummy_param = 0.0; 

    // Setup State
    State4D state = {0.0, h0, v0 * std::cos(angle_rad), v0 * std::sin(angle_rad)};
    double t = 0.0;
    double dt = 0.01; // Initial guess for time step
    
    // Statistics
    long long steps = 0;
    long long rejected_steps = 0;
    double min_dt = 1.0;
    double max_dt = 0.0;
    
    std::cout << "Running Adaptive Simulation (Tol=" << tolerance << ")...\n";
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Time(s)   X(m)      Y(m)      V(m/s)    dt(s)\n";
    std::cout << "-----------------------------------------------------\n";

    // Lambda for derivative to match signature expected by rk45_adaptive_step
    auto deriv_func = [&](const State4D& s, double time, State4D& out) {
        golf_ball_drag_deriv(s, time, out, g, dummy_param);
    };

    while (state[Y_POS] >= 0.0) {
        // 1. Attempt a step
        auto result = rk45_adaptive_step(state, t, dt, deriv_func);
        
        // 2. Calculate error magnitude (Infinity norm: max absolute error)
        double error_max = 0.0;
        for (double e : result.error) {
            error_max = std::max(error_max, std::abs(e));
        }
        
        // 3. Adaptive Controller
        if (error_max <= tolerance) {
            // -- Accept Step --
            state = result.state;
            t += dt;
            steps++;
            
            // Check dt stats
            if (dt < min_dt) min_dt = dt;
            if (dt > max_dt) max_dt = dt;

            // Output occasional status (every 0.5 sec roughly) or specific events
            // Printing every step is too much, so let's print if t crossed a 0.5s mark
            static double next_print = 0.0;
            if (t >= next_print) {
                double v = std::sqrt(state[X_VEL]*state[X_VEL] + state[Y_VEL]*state[Y_VEL]);
                std::cout << std::fixed << std::setprecision(3) 
                          << t << "     " 
                          << std::setw(8) << state[X_POS] << "  " 
                          << std::setw(8) << state[Y_POS] << "  " 
                          << std::setw(8) << v << "  " 
                          << std::scientific << std::setprecision(2) << dt << "\n";
                next_print += 0.5;
            }

            // Calculate new optimal step size
            // Formula: dt_new = dt * safety_factor * (tolerance / error)^(1/order)
            // Order is 5 for RK45 or 4 for error estimation... typically 5 or 4.
            // We use power 0.2 (1/5)
            double safety_factor = 0.9;
            double scale = safety_factor * std::pow(tolerance / (error_max + 1e-30), 0.2);
            
            // Limit growth to prevent instability
            scale = std::min(scale, 2.0);
            scale = std::max(scale, 0.1); // Don't shrink too fast either to avoid stall
            
            dt *= scale;
            
            // Cap max dt to ensure we don't skip the whole flight
            dt = std::min(dt, 0.5);

        } else {
            // -- Reject Step --
            rejected_steps++;
            
            // Shrink dt and try again
            // A simple refinement is usually half, or use the same formula as above
            dt *= 0.5;
            
            if (dt < 1e-7) {
                std::cerr << "Error: Step size underflow. Simulation failed.\n";
                break;
            }
        }
    }
    
    std::cout << "-----------------------------------------------------\n";
    std::cout << "Simulation Complete.\n";
    std::cout << "Total Distance: " << std::fixed << std::setprecision(4) << state[X_POS] << " m\n";
    std::cout << "Total Steps:    " << steps << "\n";
    std::cout << "Rejected Steps: " << rejected_steps << "\n";
    std::cout << "Min dt:         " << std::scientific << min_dt << " s\n";
    std::cout << "Max dt:         " << std::scientific << max_dt << " s\n\n";
}

int main() {
    // Run two scenarios with different tolerances to compare
    std::cout << "=== Adaptive RK45 Test: Golf Ball Trajectory ===\n\n";
    
    run_adaptive_simulation(1e-3); // Low precision (fast)
    run_adaptive_simulation(1e-6); // High precision (slower)
    
    return 0;
}
