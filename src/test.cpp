#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
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

    if (all_tests_passed) {
        std::cout << "\nAll tests passed successfully!\n";
        return 0;
    } else {
        std::cout << "\nSome tests failed.\n";
        return 1;
    }
}
