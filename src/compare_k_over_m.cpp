/**
 * @file compare_k_over_m.cpp
 * @brief Parametric study of the effect of drag properties on optimal trajectories.
 *
 * This program iterates experimentally through a range of drag coefficients
 * (k/m ratios) to observe how the optimal launch angle and maximum distance change.
 *
 * It generates a table of results showing the relationship between increasing air resistance,
 * the shallowing of the optimal launch angle, and the reduction in range.
 */
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/golden_section.hpp"
#include "../include/integrators.hpp"
#include "../include/simulation.hpp"

int main() {
    const double g = GRAVITY_EARTH;
    const double deltaT = 0.001;
    const double v0 = 100.0;
    
    // Setup integrator (RK4 is a good balance)
    SystemIntegrator<State4D> integrator = rk4_step<State4D, SystemDerivative<State4D>>;

    std::cout << "K/M Variation Analysis\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "v0 = " << v0 << " m/s\n";
    std::cout << "integrator = RK4, dt = " << deltaT << " s\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << std::setw(10) << "k/m" 
              << std::setw(25) << "Optimal Angle (deg)" 
              << std::setw(20) << "Max Distance (m)" << "\n";
    std::cout << "--------------------------------------------------------\n"; // 10 + 25 + 20 + approx spaces

    // Generate K/M values: 0.0 plus a logarithmic range from 0.001 to 0.2
    std::vector<double> k_values;
    k_values.push_back(0.0);
    
    double min_k = 0.001;
    double max_k = 0.2;
    int steps = 25; // Number of steps in the log range
    double log_min = std::log10(min_k);
    double log_max = std::log10(max_k);
    
    for (int i = 0; i <= steps; ++i) {
        double log_val = log_min + (static_cast<double>(i) / steps) * (log_max - log_min);
        k_values.push_back(std::pow(10.0, log_val));
    }

    for (double k_over_m : k_values) {
        
        // Use default derivative function (v squared drag)
        Simulation simulation(k_over_m, integrator, v0);

        auto distance_func = [&](double angle_deg) {
            return simulation.run(angle_deg, deltaT).distance;
        };

        // Standard optimization bracket [10, 80] degrees
        // Tolerance for angle finding can be reasonably tight
        double optimal_angle = golden_section_search_max(distance_func, 10.0, 80.0, 1e-4);
        double max_distance = distance_func(optimal_angle);

        std::cout << std::fixed << std::setprecision(5) 
                  << std::setw(10) << k_over_m 
                  << std::setw(25) << optimal_angle
                  << std::setw(20) << max_distance << "\n";
    }

    return 0;
}
