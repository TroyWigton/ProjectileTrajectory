/**
 * @file compare_integrators.cpp
 * @brief Analysis tool for comparing numerical integration schemes.
 *
 * This utility benchmarks different numerical integrators (Euler, Heun, RK4, RK8)
 * by measuring their performance and accuracy against a high-precision ground truth.
 *
 * It determines the required time step (dt) and number of steps for each integrator
 * to achieve a specific target distance tolerance. This helps in understanding the
 * efficiency trade-offs (computational cost vs. accuracy) for each method.
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
    const double distance_tolerance = 0.0001;
    const double k_over_m = DragRatios::PING_PONG_BALL; // Ping Pong Ball (High Drag)

    std::cout << "Project Comparison Tool\n";
    std::cout << "Drag coefficient k/m = " << k_over_m << "\n";
    std::cout << "Target precision: " << distance_tolerance << " m\n\n";

    // 1. Find optimal angle first (using RK4 as baseline)
    SystemIntegrator baseline_integrator = rk4_step<State, SystemDerivative>;
    Simulation baseline_sim(k_over_m, baseline_integrator, v0);
    auto distance_func = [&](double angle_deg) {
        return baseline_sim.run(angle_deg, 0.001).distance;
    };

    double optimal_angle = golden_section_search_max(distance_func, 10.0, 46.0, 0.001);
    std::cout << "Optimal angle used for comparison: " << optimal_angle << " degrees\n";

    // 2. Establish Ground Truth using RK8 with very small dt
    double ground_truth_dt = 1e-5;
    SystemIntegrator integrator_rk8 = rk8_step<State, SystemDerivative>;
    Simulation sim_ground_truth(k_over_m, integrator_rk8, v0);
    double ground_truth_distance = sim_ground_truth.run(optimal_angle, ground_truth_dt).distance;
    
    // Automatic precision based on tolerance (e.g. 0.0001 -> 4 decimals + 1 extra)
    int dist_precision = std::max(0, (int)ceil(-log10(distance_tolerance))) + 1;

    std::cout << "Ground Truth (RK8, dt=" << ground_truth_dt << "): " 
              << std::fixed << std::setprecision(dist_precision) << ground_truth_distance << " m\n\n";

    std::cout << "Comparison: Required Steps & dt to reach " << distance_tolerance << " m precision\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    std::cout << std::left << std::setw(15) << "Integrator" 
              << std::setw(15) << "Steps" 
              << std::setw(15) << "dt (s)" 
              << std::setw(20) << "Distance (m)" 
              << "Error (m)\n";
    std::cout << "--------------------------------------------------------------------------------\n";
    
    auto find_steps_for_precision = [&](std::string name, SystemIntegrator integ) {
         double current_dt = 1.0; 
         double current_error = 1.0e9;
         ScenarioResult result;
         
         while (true) {
             Simulation sim(k_over_m, integ, v0);
             result = sim.run(optimal_angle, current_dt);
             current_error = std::abs(result.distance - ground_truth_distance);
             
             if (current_error < distance_tolerance) break;
             
             // Refinement Strategy: Halve the time step to double the resolution
             // if the current error is too high.
             current_dt /= 2.0;

             // Safety Break: Prevent infinite loops or excessive computation
             // if the integrator fails to converge within a reasonable limit.
             if (current_dt < 1e-7) break;
         }


        std::cout << std::left << std::setw(15) << name 
                  << std::setw(15) << result.step_count 
                  << std::setw(15) << std::defaultfloat << std::setprecision(6) << result.deltaT 
                  << std::setw(20) << std::fixed << std::setprecision(dist_precision) << result.distance 
                  << std::scientific << std::setprecision(2) << current_error << "\n";
    };

    find_steps_for_precision("Euler", euler_step<State, SystemDerivative>);
    find_steps_for_precision("Heun", heun_step<State, SystemDerivative>);
    find_steps_for_precision("RK4", rk4_step<State, SystemDerivative>);
    find_steps_for_precision("RK8", rk8_step<State, SystemDerivative>);

    return 0;
}
