/**
 * @file compare_integrators.cpp
 * @brief Analysis tool for comparing numerical integration schemes.
 *
 * This utility benchmarks different numerical integrators (Euler, Heun, RK4, RK45, RK8)
 * by measuring their performance and accuracy against a high-precision ground truth.
 *
 * It determines the required time step (dt) and number of steps for each integrator
 * to achieve a specific target precision (Position Error Norm at a fixed time T).
 * This helps in understanding the efficiency trade-offs (computational cost vs. accuracy)
 * for each method.
 */
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/golden_section.hpp"
#include "../include/integrators.hpp"
#include "../include/simulation.hpp"

int main() {
    const double g = GRAVITY_EARTH;
    const double v0 = 100;
    const double h0 = 0.0;
    const double distance_tolerance = 0.001; // 1.0 mm (Common Engineering Tolerance)
    const double k_over_m = DragRatios::PING_PONG_BALL; // Higher drag = "harder" problem
    const double t_end = 3.0; // Shorter duration to allow Euler to verify easier

    std::cout << "Numerical Integrator Benchmark\n";
    std::cout << "Drag/Mass Ratio (k/m) = " << k_over_m << "\n";
    std::cout << "Initial Velocity (v0) = " << v0 << " m/s\n";
    std::cout << "Target Precision: " << distance_tolerance << " m (Position Error Norm at fixed time T=" << t_end << "s)\n\n";

    // 1. Establish Ground Truth using RK8 with very small dt
    // We use a lambda to run a fixed-time simulation
    auto run_fixed_time = [&](SystemIntegrator<State4D> integrator, double dt) -> State4D {
        using namespace StateIndex4D;
        double angle_rad = 45.0 * M_PI / 180.0; // Fixed angle
        State4D state = {0.0, h0, v0 * std::cos(angle_rad), v0 * std::sin(angle_rad)};
        double t = 0.0;
        
        // Lambda to bind parameters
        auto bound_deriv = [&](const State4D& s, double time, State4D& out_d) {
             drag_deriv_v_squared(s, time, out_d, g, k_over_m); 
        };

        long long steps = std::round(t_end / dt);
        // Recalculate exact dt to hit t_end precisely
        double actual_dt = t_end / steps; 

        for(long long i=0; i<steps; ++i) {
            state = integrator(state, t, actual_dt, bound_deriv);
            t += actual_dt;
        }
        return state;
    };

    const double ground_truth_dt = 1e-5;
    State4D ground_truth_state = run_fixed_time(rk8_step<State4D, SystemDerivative<State4D>>, ground_truth_dt);
    
    std::cout << "Reference Truth (RK8, dt=" << ground_truth_dt << "s) at t=" << t_end << "s: (" 
              << ground_truth_state[0] << ", " << ground_truth_state[1] << ")\n\n";

    std::cout << "Methodology: Iteratively reducing time step (dt) until error < target precision.\n";
    std::cout << "Performance Comparison: Steps & dt required to meet target precision\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << std::left << std::setw(15) << "Integrator" 
              << std::setw(15) << "Steps" 
              << std::setw(15) << "dt (s)" 
              << std::setw(20) << "Final Error (m)" 
              << "\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    auto find_steps_for_precision = [&](std::string name, SystemIntegrator<State4D> integ) {
         double current_dt = 1.0; 
         double current_error = 1.0e9;
         State4D result;
         long long step_count = 0;
         
         while (true) {
             result = run_fixed_time(integ, current_dt);
             step_count = std::round(t_end / current_dt);

             // Calculate Euclidean distance in position only (x,y)
             double dx = result[0] - ground_truth_state[0];
             double dy = result[1] - ground_truth_state[1];
             current_error = std::sqrt(dx*dx + dy*dy);
             
             if (current_error < distance_tolerance) break;
             
             // Refinement Strategy: Reduce time step gradually to find limit
             current_dt *= 0.99;

             // Safety Break: Prevent infinite loops or excessive computation
             if (step_count > 20000000) break; 
         }

        // Recalculate dt for display based on final step count to match run_fixed_time logic
        double final_dt = t_end / step_count;

        std::cout << std::left << std::setw(15) << name 
                  << std::setw(15) << step_count 
                  << std::setw(15) << std::defaultfloat << std::setprecision(6) << final_dt 
                  << std::scientific << std::setprecision(2) << current_error << "\n";
    };

    find_steps_for_precision("Euler", euler_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("Heun", heun_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("RK4", rk4_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("RK45", rk45_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("RK8", rk8_step<State4D, SystemDerivative<State4D>>);

    return 0;
}
