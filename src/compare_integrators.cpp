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
#include <sstream>
#include <chrono>
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
    const double v0 = 200.0;
    const double h0 = 0.0;
    const double launch_angle_deg = 45.0;
    const double distance_tolerance = 0.001; // 1.0 mm (Common Engineering Tolerance)
    const double k_over_m = DragRatios::PING_PONG_BALL;
    const double t_end = 20.0; // Longer duration => more steps => runtimes rise above clock jitter

    std::cout << "Numerical Integrator Benchmark\n";
    std::cout << "Simulates a 2D projectile (gravity + v^2 air drag) for a fixed duration T,\n";
    std::cout << "then compares each integrator's final (x, y) position against an RK8 ground\n";
    std::cout << "truth to measure accuracy vs. computational cost.\n\n";
    std::cout << "Launch Angle:         " << launch_angle_deg << " deg\n";
    std::cout << "Initial Velocity (v0): " << v0 << " m/s\n";
    std::cout << "Drag/Mass Ratio (k/m): " << k_over_m << "\n";
    std::cout << "Simulation Duration:  " << t_end << " s\n";
    std::cout << "Target Precision:     " << distance_tolerance << " m (Position Error Norm at t=T)\n\n";

    // 1. Establish Ground Truth using RK8 with very small dt
    // We use a lambda to run a fixed-time simulation
    auto run_fixed_time = [&](SystemIntegrator<State4D> integrator, double dt) -> State4D {
        using namespace StateIndex4D;
        double angle_rad = launch_angle_deg * M_PI / 180.0;
        State4D state = {0.0, h0, v0 * std::cos(angle_rad), v0 * std::sin(angle_rad)};
        double t = 0.0;
        
        // Lambda to bind parameters via a ProjectileContext
        ProjectileContext ctx{g, k_over_m};
        auto bound_deriv = [&ctx](const State4D& s, double time, State4D& out_d) {
             drag_deriv_v_squared(s, time, out_d, ctx);
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
    
    std::cout << "Reference Truth (RK8, dt=" << ground_truth_dt << "s) at t=" << t_end << "s: (x="
              << ground_truth_state[0] << " m, y=" << ground_truth_state[1] << " m)\n\n";

    std::cout << "Methodology: Iteratively reducing time step (dt) until position error < target precision.\n";
    std::cout << "             Runtime is measured by re-running once at the converged dt.\n";
    std::cout << "Performance Comparison: Steps, dt, error, and wall-clock runtime to meet target precision\n";
    std::cout << "--------------------------------------------------------------------------------------------\n";
    std::cout << std::left << std::setw(15) << "Integrator"
              << std::setw(15) << "Steps"
              << std::setw(15) << "dt (s)"
              << std::setw(20) << "Final Error (m)"
              << "Runtime (ms)\n";
    std::cout << "--------------------------------------------------------------------------------------------\n";

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

        // Timed re-run at the converged dt to measure wall-clock cost of step_count steps.
        // The volatile sink forces the optimizer to keep the simulation result, preventing
        // dead-code elimination of the timed call under -O3.
        auto t0 = std::chrono::high_resolution_clock::now();
        State4D timed_result = run_fixed_time(integ, final_dt);
        auto t1 = std::chrono::high_resolution_clock::now();
        volatile double sink = timed_result[0] + timed_result[1];
        (void)sink;
        double runtime_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // Stream error into a string so setw() controls its column width cleanly.
        std::ostringstream err_str;
        err_str << std::scientific << std::setprecision(2) << current_error;

        std::cout << std::left << std::setw(15) << name
                  << std::setw(15) << step_count
                  << std::setw(15) << std::defaultfloat << std::setprecision(6) << final_dt
                  << std::setw(20) << err_str.str()
                  << std::fixed << std::setprecision(3) << runtime_ms << "\n";
    };

    find_steps_for_precision("Euler", euler_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("Heun", heun_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("RK4", rk4_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("RK45", rk45_step<State4D, SystemDerivative<State4D>>);
    find_steps_for_precision("RK8", rk8_step<State4D, SystemDerivative<State4D>>);

    return 0;
}
