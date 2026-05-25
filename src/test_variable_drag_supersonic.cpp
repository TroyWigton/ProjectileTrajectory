/**
 * @file test_variable_drag_supersonic.cpp
 * @brief Compare constant-Cd vs variable-Cd (Mach-dependent) drag on a supersonic projectile.
 *
 * Evaluates the impact of the transonic drag rise and supersonic fade on a high-velocity
 * projectile (.223 Remington, ~Mach 2.84 at the muzzle). Trajectories compared:
 *   Model 0: Vacuum (no drag) — theoretical baseline
 *   Model 1: v^2 drag with constant Cd
 *   Model 2: v^2 drag with variable Cd (transonic rise with supersonic fade)
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

    std::cout << "Supersonic Variable Drag Model Evaluation (.223 Remington)\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "Initial Velocity: " << v0 << " m/s (Mach " << std::fixed << std::setprecision(4) << v0/343.0 << ")\n";
    std::cout << "Base k/m: " << k_over_m << "\n";
    std::cout << "Time step: " << deltaT << " s\n";
    std::cout << "--------------------------------------------------------\n\n";

    // Setup Integrator (RK4)
    SystemIntegrator<State4D> integrator = rk4_step<State4D, SystemDerivative<State4D>>;

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

    // 1. v^2 drag with constant Cd
    Simulation sim_standard(k_over_m, integrator, v0, drag_deriv_v_squared, g, h0);

    auto dist_func_std = [&](double angle) {
        return sim_standard.run(angle, deltaT).distance;
    };

    double opt_angle_std = golden_section_search_max(dist_func_std, 10.0, 60.0, 1e-4);
    double max_dist_std = dist_func_std(opt_angle_std);

    std::cout << "Model 1: v^2 drag, constant Cd\n";
    std::cout << "  Optimal Angle: " << std::fixed << std::setprecision(4) << opt_angle_std << " degrees\n";
    std::cout << "  Max Distance:  " << std::fixed << std::setprecision(2) << max_dist_std << " m\n\n";


    // 2. v^2 drag with variable Cd (Mach-dependent: transonic rise with supersonic fade)
    Simulation sim_variable(k_over_m, integrator, v0, variable_drag_deriv, g, h0);

    auto dist_func_var = [&](double angle) {
        return sim_variable.run(angle, deltaT).distance;
    };

    double opt_angle_var = golden_section_search_max(dist_func_var, 10.0, 60.0, 1e-4);
    double max_dist_var = dist_func_var(opt_angle_var);

    std::cout << "Model 2: v^2 drag, variable Cd (transonic rise with supersonic fade)\n";
    std::cout << "  Optimal Angle: " << std::fixed << std::setprecision(4) << opt_angle_var << " degrees\n";
    std::cout << "  Max Distance:  " << std::fixed << std::setprecision(2) << max_dist_var << " m\n\n";

    // Compare Model 1 (constant Cd) vs Model 2 (variable Cd) to isolate the
    // transonic/supersonic Cd variation effect.
    double dist_diff = max_dist_std - max_dist_var;
    double percent_diff = (dist_diff / max_dist_std) * 100.0;

    std::cout << "Impact Analysis (Model 2 variable Cd vs Model 1 constant Cd):\n";
    std::cout << "  Range Reduction: " << dist_diff << " m (" << percent_diff << "% loss vs constant Cd)\n";
    std::cout << "  Angle Shift:     " << (opt_angle_var - opt_angle_std) << " degrees (variable Cd - constant Cd)\n";

    return 0;
}
