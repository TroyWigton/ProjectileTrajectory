/**
 * @file test_variable_drag.cpp
 * @brief Test suite for comparing constant vs variable (Mach-dependent) drag models.
 *
 * This test evaluates the impact of transonic/supersonic drag rise on a high-velocity
 * projectile (.223 Remington). It compares the trajectories predicted by:
 * 1. Standard Velocity-Squared Drag (constant k/m)
 * 2. Variable Drag (k/m changes with Mach number)
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/golden_section.hpp"
#include "../include/integrators.hpp"
#include "../include/simulation.hpp"

int main() {
    const double g = GRAVITY_EARTH;
    const double deltaT = 0.001;
    // High velocity projectile: .223 Remington
    const double v0 = 975.0; // m/s
    const double h0 = 0.0;
    
    // Base k/m derived from Cd=0.23 for a 55gr bullet
    const double k_over_m = DragRatios::BULLET_223;

    std::cout << "Variable Drag Model Evaluation (.223 Remington)\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "Initial Velocity: " << v0 << " m/s (Mach " << std::fixed << std::setprecision(4) << v0/343.0 << ")\n";
    std::cout << "Base k/m: " << k_over_m << "\n";
    std::cout << "Time step: " << deltaT << " s\n";
    std::cout << "--------------------------------------------------------\n\n";

    // Setup Integrator (RK4)
    SystemIntegrator integrator = rk4_step<State, SystemDerivative>;

    // 0. Vacuum Model (No Drag)
    // Pass 0.0 for drag coefficient, though no_drag_deriv ignores it anyway
    Simulation sim_vacuum(0.0, integrator, v0, no_drag_deriv, g, h0);
    
    // For vacuum, we know optimal is 45 degrees, but let's calculate it to be consistent
    auto dist_func_vac = [&](double angle) {
        return sim_vacuum.run(angle, deltaT).distance;
    };

    double opt_angle_vac = golden_section_search_max(dist_func_vac, 10.0, 60.0, 1e-4);
    double max_dist_vac = dist_func_vac(opt_angle_vac);

    std::cout << "Model 0: Vacuum (Theoretical Baseline)\n";
    std::cout << "  Optimal Angle: " << std::fixed << std::setprecision(4) << opt_angle_vac << " degrees\n";
    std::cout << "  Max Distance:  " << std::fixed << std::setprecision(2) << max_dist_vac << " m\n\n";

    // 1. Standard Drag Model (Constant k/m)
    Simulation sim_standard(k_over_m, integrator, v0, drag_deriv_v_squared, g, h0);
    
    auto dist_func_std = [&](double angle) {
        return sim_standard.run(angle, deltaT).distance;
    };

    double opt_angle_std = golden_section_search_max(dist_func_std, 10.0, 60.0, 1e-4);
    double max_dist_std = dist_func_std(opt_angle_std);
    
    std::cout << "Model 1: Standard v^2 Drag (Constant Coeff)\n";
    std::cout << "  Optimal Angle: " << std::fixed << std::setprecision(4) << opt_angle_std << " degrees\n";
    std::cout << "  Max Distance:  " << std::fixed << std::setprecision(2) << max_dist_std << " m\n\n";


    // 2. Variable Drag Model (Mach Dependent)
    Simulation sim_variable(k_over_m, integrator, v0, variable_drag_deriv, g, h0);
    
    auto dist_func_var = [&](double angle) {
        return sim_variable.run(angle, deltaT).distance;
    };

    double opt_angle_var = golden_section_search_max(dist_func_var, 10.0, 60.0, 1e-4);
    double max_dist_var = dist_func_var(opt_angle_var);

    std::cout << "Model 2: Variable Drag (Transonic Rise)\n";
    std::cout << "  Optimal Angle: " << std::fixed << std::setprecision(4) << opt_angle_var << " degrees\n";
    std::cout << "  Max Distance:  " << std::fixed << std::setprecision(2) << max_dist_var << " m\n\n";

    // Comparison
    double dist_diff = max_dist_std - max_dist_var;
    double percent_diff = (dist_diff / max_dist_std) * 100.0;
    
    std::cout << "Impact Analysis:\n";
    std::cout << "  Range Reduction: " << dist_diff << " m (" << percent_diff << "% loss)\n";
    std::cout << "  Angle Shift:     " << (opt_angle_var - opt_angle_std) << " degrees\n";

    return 0;
}
