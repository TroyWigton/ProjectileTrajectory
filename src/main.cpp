#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <functional>
#include <string>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"

//#define DEBUG // define before local includes to enable debug mode in those files
#include "../include/golden_section.hpp"
#include "../include/integrators.hpp"
#include "../include/simulation.hpp"

int main() {
    const double g = 9.81;
    const double deltaT = 0.001;
    const double v0 = 100;
    const double h0 = 0.0;
    const double distance_tolerance = 0.01; // Maximum distance loss (m) tolerable at angle precision boundary

    double k_over_m; // 
    // Vacuum 0.0, GolfBall .0025, PingPongBall .01
    // k_over_m = 0.0;   // No Drag scenario (validation testing)
    k_over_m = 0.0025;
    std::cout << "Using drag coefficient k/m = " << k_over_m << "\n";
    std::cout << "Target distance precision: " << distance_tolerance << " m\n";

    // Cast the template function to the specific SystemIntegrator type
    // to allow implicit conversion in Simulation constructor
    SystemIntegrator integrator = rk4_step<State, SystemDerivative>;

    // Use default drag_deriv_v_squared
    Simulation simulation(k_over_m, integrator, v0);

    auto distance_func = [&](double angle_deg) {
        return simulation.run(angle_deg, deltaT).distance;
    };

    // Create a wide initial bracket that definitely contains the maximum
    double a = 10.0;
    double b = 46.0;
    std::cout << std::fixed << std::setprecision(6);

    // Start with a coarse angle tolerance to find the optimal angle quickly
    double coarse_angle_tol = 0.01;
    std::cout << "\nPhase 1: Finding approximate optimal angle (coarse search):\n";
    double optimal_angle = golden_section_search_max(distance_func, a, b, coarse_angle_tol);
    double max_distance = distance_func(optimal_angle);
    std::cout << "Approximate optimal angle: " << optimal_angle << " degrees\n";
    std::cout << "Maximum distance: " << max_distance << " m\n";

    // Now refine to find the angle precision needed for the distance tolerance
    std::cout << "\nPhase 2: Finding angle precision for distance tolerance " 
              << distance_tolerance << " m:\n";
    
    // Binary search to find the angle offset that produces the distance tolerance
    double angle_offset = 0.1;  // Start with 0.1 degree offset
    double low = 0.0;
    double high = 5.0;  // Maximum reasonable angle deviation
    
    while (high - low > 1e-9) {  // Very fine precision for angle determination
        angle_offset = (low + high) / 2.0;
        
        double dist_plus = distance_func(optimal_angle + angle_offset);
        double dist_minus = distance_func(optimal_angle - angle_offset);
        
        // Check if either side falls outside the tolerance
        double max_deviation = std::max(
            std::abs(max_distance - dist_plus),
            std::abs(max_distance - dist_minus)
        );
        
        if (max_deviation > distance_tolerance) {
            high = angle_offset;  // Angle offset is too large
        } else {
            low = angle_offset;   // Angle offset might be acceptable
        }
    }
    
    double required_angle_precision = angle_offset;
    
    // Now refine the optimal angle to this precision
    std::cout << "\nPhase 3: Refining optimal angle to required precision:\n";
    a = optimal_angle - 2.0 * required_angle_precision;
    b = optimal_angle + 2.0 * required_angle_precision;
    optimal_angle = golden_section_search_max(distance_func, a, b, required_angle_precision / 10.0);
    max_distance = distance_func(optimal_angle);
    
    // Calculate appropriate precision for displaying angle precision
    // Find the order of magnitude of the angle precision
    int angle_precision_digits = 1;
    if (required_angle_precision > 0) {
        angle_precision_digits = std::max(1, static_cast<int>(-std::log10(required_angle_precision)) + 1);
    }
    
    std::cout << "\n=== RESULTS ===\n";
    std::cout << "Optimal launch angle: " << optimal_angle << " degrees\n";
    std::cout << "Maximum distance: " << max_distance << " m\n";
    std::cout << std::setprecision(angle_precision_digits);
    std::cout << "Required angle precision: ±" << required_angle_precision << " degrees\n";
    std::cout << std::setprecision(6);  // Reset for other outputs
    std::cout << "  (to maintain distance within ±" << distance_tolerance << " m)\n";
    
    // Final verification
    std::cout << "\nVerification at angle tolerance boundaries:\n";
    double dist_at_optimal = distance_func(optimal_angle);
    double dist_plus = distance_func(optimal_angle + required_angle_precision);
    double dist_minus = distance_func(optimal_angle - required_angle_precision);
    
    std::cout << "  Angle: " << (optimal_angle - required_angle_precision) 
              << "° → Distance: " << dist_minus 
              << " m (Δ = " << (max_distance - dist_minus) << " m)\n";
    std::cout << "  Angle: " << optimal_angle 
              << "° → Distance: " << dist_at_optimal << " m (optimal)\n";
    std::cout << "  Angle: " << (optimal_angle + required_angle_precision) 
              << "° → Distance: " << dist_plus 
              << " m (Δ = " << (max_distance - dist_plus) << " m)\n";

    return 0;
}