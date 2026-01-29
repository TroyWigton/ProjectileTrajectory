/**
 * @file test.cpp
 * @brief Unit tests and validation suite for the Projectile Simulation library.
 *
 * This file contains a comprehensive test suite to verify the correctness and
 * accuracy of the numerical integrators and trajectory simulations.
 *
 * Tests included:
 * 1. Vacuum Trajectory Optimization: Verifies that the optimal angle is 45 degrees
 *    when no drag is present (theoretical truth).
 * 2. Integrator Consistency: Compares the results of Euler, Heun, and RK4 methods
 *    against a high-precision Runge-Kutta 8 (RK8) baseline under drag conditions.
 *    Ensures that less accurate methods drift within expected bounds relative to
 *    simpler methods.
 * 3. Vacuum Trajectory with initial height: Verifies that the simulated optimal angle matches
 *    the closed-form theoretical prediction when launching from a non-zero initial height.
 */
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/golden_section.hpp"
#include "../include/integrators.hpp"
#include "../include/simulation.hpp"

// Define a struct to hold integrator information
struct IntegratorInfo {
    std::string name;
    SystemIntegrator func;
};

// Define a struct to hold derivative function information
struct DerivativeInfo {
    std::string name;
    DerivativeFuncPtr func;
};

int main() {
    const double g = 9.81;
    const double deltaT = 0.001;
    const double v0 = 100.0;
    const double h0 = 0.0;
    const double k_over_m = 0.0; // Vacuum

    std::vector<IntegratorInfo> integrators = {
        {"Euler", euler_step<State, SystemDerivative>},
        {"Heun", heun_step<State, SystemDerivative>},
        {"RK4", rk4_step<State, SystemDerivative>},
        {"RK8", rk8_step<State, SystemDerivative>}
    };

    std::vector<DerivativeInfo> derivatives = {
        {"V Squared Drag", drag_deriv_v_squared},
        {"Linear Drag", drag_deriv_linear},
        {"No Drag", no_drag_deriv}
    };

    bool all_tests_passed = true;
    const double expected_angle = 45.0;
    const double tolerance = 0.05; // Tolerance for angle comparison

    std::cout << "Starting Test Suite: Vacuum Trajectory Optimization\n";
    std::cout << "Target Angle: " << expected_angle << " degrees\n";
    std::cout << "Tolerance: " << tolerance << " degrees\n";
    std::cout << "k/m: " << k_over_m << "\n\n";

    for (const auto& integ : integrators) {
        for (const auto& deriv : derivatives) {
            std::cout << "Testing Integrator: " << integ.name << ", Derivative: " << deriv.name << "... ";

            Simulation simulation(k_over_m, integ.func, v0, deriv.func, g, h0);

            auto distance_func = [&](double angle_deg) {
                return simulation.run(angle_deg, deltaT).distance;
            };

            // Bracket for search
            double a = 10.0;
            double b = 80.0;
            double angle_tol = 1e-4;

            double optimal_angle = golden_section_search_max(distance_func, a, b, angle_tol);

            if (std::abs(optimal_angle - expected_angle) < tolerance) {
                std::cout << "PASSED (Optimal Angle: " << optimal_angle << ")\n";
            } else {
                std::cout << "FAILED (Optimal Angle: " << optimal_angle << ")\n";
                all_tests_passed = false;
            }
        }
    }

    std::cout << "\n--------------------------------------------------\n";
    std::cout << "Starting Test Suite: Integrator Consistency (Non-zero Drag)\n";
    double k_over_m_drag = 0.0057; 
    std::cout << "k/m: " << k_over_m_drag << "\n";
    
    // We expect RK8 to be the most accurate, valid baseline.
    // We will compare others against RK8.
    
    for (const auto& deriv : derivatives) {
        if (deriv.name == "No Drag") continue; // Skip no drag for this test
        
        std::cout << "\nEvaluating Derivative Function: " << deriv.name << "\n";
        
        // Store results to compare
        struct Result {
            double angle;
            double distance;
        };
        std::map<std::string, Result> results;
        
        // Run for each integrator
        for (const auto& integ : integrators) {
             Simulation simulation(k_over_m_drag, integ.func, v0, deriv.func, g, h0);
             auto distance_func = [&](double angle_deg) {
                return simulation.run(angle_deg, deltaT).distance;
            };
            
            double optimal_angle = golden_section_search_max(distance_func, 10.0, 80.0, 1e-4);
            double max_dist = distance_func(optimal_angle);
            
            results[integ.name] = {optimal_angle, max_dist};
            std::cout << "  " << std::setw(10) << integ.name 
                      << ": Angle = " << optimal_angle 
                      << ", Distance = " << max_dist << "\n";
        }
        
        // Verification: Compare all against RK8 baseline
        if (results.find("RK8") != results.end()) {
            const auto& rk8_result = results["RK8"];
            
            // Define tolerances for each method (Angle, Distance)
            // Euler (1st order) needs looser tolerances than RK4 (4th order)
            std::map<std::string, std::pair<double, double>> tolerances = {
                {"Euler", {0.1, 0.5}},    // Coarser approximation
                {"Heun",  {0.01, 0.1}},   // 2nd order accuracy
                {"RK4",   {0.001, 0.05}}  // High accuracy
            };

            for (const auto& entry : results) {
                const std::string& name = entry.first;
                if (name == "RK8") continue; // Don't compare with self

                double angle_diff = std::abs(entry.second.angle - rk8_result.angle);
                double dist_diff = std::abs(entry.second.distance - rk8_result.distance);
                
                double angle_tol = tolerances[name].first;
                double dist_tol = tolerances[name].second;

                if (angle_diff < angle_tol && dist_diff < dist_tol) {
                    std::cout << "  -> Consistency Check (" << name << " vs RK8): PASSED\n";
                    std::cout << "     (Diffs - Angle: " << angle_diff << ", Dist: " << dist_diff << ")\n";
                } else {
                    std::cout << "  -> Consistency Check (" << name << " vs RK8): FAILED\n";
                    std::cout << "     (Diffs - Angle: " << angle_diff << ", Dist: " << dist_diff 
                              << " | Tolerance: " << angle_tol << ", " << dist_tol << ")\n";
                    all_tests_passed = false;
                }
            }
        }
    }

    std::cout << "\n--------------------------------------------------\n";
    std::cout << "Starting Test Suite: Vacuum Trajectory with Height (No Drag)\n";
    
    double h0_test = 100.0; // Launch height
    std::cout << "Launch Height h0: " << h0_test << " m\n";
    std::cout << "Initial Velocity v0: " << v0 << " m/s\n";
    std::cout << "Gravity g: " << g << " m/s^2\n";
    
    // Theoretical Calculation
    // Theta = arctan( v0 / sqrt(v0^2 + 2*g*h0) )
    double theta_rad = std::atan(v0 / std::sqrt(v0 * v0 + 2 * g * h0_test));
    double theoretical_angle = theta_rad * 180.0 / M_PI;
    
    std::cout << "Theoretical Optimal Angle: " << theoretical_angle << " degrees\n";
    
    // Simulation
    // Use RK4 and no_drag_deriv
    SystemIntegrator integ_rk4 = rk4_step<State, SystemDerivative>;
    // k_over_m = 0.0 for vacuum
    Simulation simulation_h(0.0, integ_rk4, v0, no_drag_deriv, g, h0_test);
    
    auto distance_func_h = [&](double angle_deg) {
        return simulation_h.run(angle_deg, deltaT).distance;
    };
    
    // Search around the theoretical value
    double search_range = 10.0;
    double experimental_angle_h = golden_section_search_max(distance_func_h, 
                                                           theoretical_angle - search_range, 
                                                           theoretical_angle + search_range, 
                                                           1e-4);
    
    std::cout << "Experimental Optimal Angle: " << experimental_angle_h << " degrees\n";
    
    if (std::abs(experimental_angle_h - theoretical_angle) < tolerance) {
        std::cout << "PASSED (Difference: " << std::abs(experimental_angle_h - theoretical_angle) << ")\n";
    } else {
        std::cout << "FAILED (Difference: " << std::abs(experimental_angle_h - theoretical_angle) << ")\n";
        all_tests_passed = false;
    }

    // Distance Verification
    // R_max = (v0 / g) * sqrt(v0^2 + 2*g*h0)
    double theoretical_max_dist = (v0 / g) * std::sqrt(v0 * v0 + 2 * g * h0_test);
    double experimental_max_dist = distance_func_h(experimental_angle_h);
    double dist_tolerance = 0.05; // 5cm tolerance

    std::cout << "Theoretical Max Distance: " << theoretical_max_dist << " m\n";
    std::cout << "Experimental Max Distance: " << experimental_max_dist << " m\n";

    if (std::abs(experimental_max_dist - theoretical_max_dist) < dist_tolerance) {
        std::cout << "Distance Check: PASSED (Difference: " << std::abs(experimental_max_dist - theoretical_max_dist) << ")\n";
    } else {
        std::cout << "Distance Check: FAILED (Difference: " << std::abs(experimental_max_dist - theoretical_max_dist) << ")\n";
        all_tests_passed = false;
    }

    if (all_tests_passed) {
        std::cout << "\nAll tests passed successfully!\n";
        return 0;
    } else {
        std::cout << "\nSome tests failed.\n";
        return 1;
    }
}
