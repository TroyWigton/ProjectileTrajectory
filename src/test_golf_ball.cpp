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

    std::cout << std::left << std::setw(26) << label
              << " Optimal Angle: " << std::fixed << std::setprecision(4) << optimal_angle << " deg"
              << " | Max Distance: " << std::setprecision(3) << max_dist << " m\n";
}

int main() {
    std::cout << "Golf Ball Trajectory Optimization (v0=70 m/s ~ 156 mph ball speed)\n";
    std::cout << "------------------------------------------------------------------------\n";

    double v0 = 70.0; // Typical driver ball speed

    // Cd values used by cases 2 and 3. k/m = F * Cd, where F = (rho * A) / (2 * m).
    const double cd_laminar    = GolfBallPhysics::CD_LAMINAR;
    const double cd_turbulent  = GolfBallPhysics::CD_TURBULENT;
    const double cd_mean       = (cd_laminar + cd_turbulent) / 2.0;
    const double v_trans_start = GolfBallPhysics::V_TRANSITION_START;
    const double v_trans_end   = GolfBallPhysics::V_TRANSITION_END;
    const double km_std        = GolfBallPhysics::FACTOR_F * cd_mean;

    std::cout << std::fixed;
    std::cout << "Cd model parameters:\n";
    std::cout << "  Constant Cd: " << std::setprecision(4) << cd_mean
              << "  (mean of laminar=" << std::setprecision(3) << cd_laminar
              << " and turbulent=" << cd_turbulent << ")\n";
    std::cout << "  Variable Cd: " << cd_laminar << " (laminar, v < "
              << std::setprecision(1) << v_trans_start << " m/s) -> "
              << std::setprecision(3) << cd_turbulent
              << " (turbulent, v >= " << std::setprecision(1) << v_trans_end << " m/s)\n";
    std::cout << "               with linear transition between "
              << v_trans_start << " and " << v_trans_end << " m/s\n";
    std::cout << "------------------------------------------------------------------------\n";

    // 1. No drag (vacuum baseline)
    run_optimization_case("1. Vacuum", v0, 0.0, no_drag_deriv);

    // 2. v^2 drag with a constant Cd (mean of laminar and turbulent values)
    run_optimization_case("2. v^2 drag, constant Cd", v0, km_std, drag_deriv_v_squared);

    // 3. v^2 drag with velocity-dependent Cd.
    //    k_over_m argument is ignored; the model uses GolfBallPhysics constants directly.
    run_optimization_case("3. v^2 drag, variable Cd", v0, 0.0, golf_ball_drag_deriv);

    return 0;
}
