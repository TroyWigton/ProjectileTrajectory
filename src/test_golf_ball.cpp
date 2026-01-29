#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/integrators.hpp"
#include "../include/simulation.hpp"
#include "../include/golden_section.hpp"

// Optimize angle for a given configuration and print results
void run_optimization_case(std::string label, double v0, double k_over_m, DerivativeFuncPtr deriv_func) {
    const double g = GRAVITY_EARTH;
    const double dt = 0.001; 
    const double h0 = 0.0;

    // Use standard RK4
    SystemIntegrator<State4D> integrator = rk4_step<State4D, SystemDerivative<State4D>>;

    // Create Simulation instance
    Simulation sim(k_over_m, integrator, v0, deriv_func, g, h0);

    // Lambda for optimization
    auto distance_func = [&](double angle_deg) {
        return sim.run(angle_deg, dt).distance;
    };

    // Find optimal angle using Golden Section Search
    // Bracket: 5 to 55 degrees should cover all reasonable golf trajectories
    double optimal_angle = golden_section_search_max(distance_func, 5.0, 55.0, 1e-4);
    double max_dist = distance_func(optimal_angle);

    std::cout << std::left << std::setw(20) << label 
              << " Optimal Angle: " << std::fixed << std::setprecision(4) << optimal_angle << " deg"
              << " | Max Distance: " << std::setprecision(3) << max_dist << " m\n";
}

int main() {
    std::cout << "Golf Ball Trajectory Optimization (v0=70 m/s ~ 156 mph ball speed)\n";
    std::cout << "------------------------------------------------------------------------\n";
    
    double v0 = 70.0; // Typical driver ball speed
    
    // 1. Vacuum
    run_optimization_case("1. Vacuum", v0, 0.0, no_drag_deriv);

    // 2. Constant Drag (Baseline using Cd_final / Turbulent)
    // k/m = F * Cd_final
    double km_std = GolfBallPhysics::FACTOR_F * GolfBallPhysics::CD_LAMINAR;
    run_optimization_case("2. Const Low Drag", v0, km_std, drag_deriv_v_squared);

    // 3. New Golf Ball Model (Variable Drag)
    run_optimization_case("3. Variable Drag", v0, 0.0, golf_ball_drag_deriv); // params ignored, uses constants

    return 0;
}
