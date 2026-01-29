#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/integrators.hpp"

// Define a 6D derivative function for 3D motion with drag
void projectile_3d_deriv(const State6D& s, double t, State6D& deriv, double g, double k_over_m) {
    using namespace StateIndex6D;

    // Position derivatives (velocities)
    deriv[X_POS] = s[X_VEL];
    deriv[Y_POS] = s[Y_VEL];
    deriv[Z_POS] = s[Z_VEL];

    // Calculate velocity magnitude
    double v_sq = s[X_VEL]*s[X_VEL] + s[Y_VEL]*s[Y_VEL] + s[Z_VEL]*s[Z_VEL];
    double v = std::sqrt(v_sq);

    // Forces (Drag + Gravity)
    // Drag acts opposite to velocity vector
    // a_drag_x = -(k/m) * v * vx
    deriv[X_VEL] = -k_over_m * v * s[X_VEL];
    deriv[Y_VEL] = -k_over_m * v * s[Y_VEL] - g; // Gravity acts on Y
    deriv[Z_VEL] = -k_over_m * v * s[Z_VEL];
}

int main() {
    using namespace StateIndex6D;

    double g = 9.81;
    double k_over_m = 0.005; // Low drag
    double v0 = 100.0;
    double angle_launch = 45.0 * M_PI / 180.0;
    double angle_side = 10.0 * M_PI / 180.0; // Launch slightly sideways

    // Initial state: x, y, z, vx, vy, vz
    State6D state;
    state[X_POS] = 0.0;
    state[Y_POS] = 0.0;
    state[Z_POS] = 0.0;
    state[X_VEL] = v0 * std::cos(angle_launch) * std::cos(angle_side);
    state[Y_VEL] = v0 * std::sin(angle_launch);
    state[Z_VEL] = v0 * std::cos(angle_launch) * std::sin(angle_side);

    // Setup integrator
    // We can use the generic rk4_step directly
    SystemIntegrator<State6D> integrator = rk4_step<State6D, std::function<void(const State6D&, double, State6D&)>>;

    // Bind parameters to derivative function
    auto deriv_func = [&](const State6D& s, double t, State6D& d) {
        projectile_3d_deriv(s, t, d, g, k_over_m);
    };

    double dt = 0.01;
    double t = 0.0;
    
    std::cout << "Running 3D Projectile Simulation (6-DOF State)...\n";
    std::cout << "Drag (k/m): " << k_over_m << "\n";
    std::cout << std::fixed << std::setprecision(4);

    while (state[Y_POS] >= 0.0) {
        state = integrator(state, t, dt, deriv_func);
        t += dt;

        // Print every 1 second
        if (std::abs(std::fmod(t, 1.0)) < dt/2.0) {
            std::cout << "t=" << t 
                      << " x=" << state[X_POS] 
                      << " y=" << state[Y_POS] 
                      << " z=" << state[Z_POS] 
                      << " v=" << std::sqrt(state[X_VEL]*state[X_VEL] + state[Y_VEL]*state[Y_VEL] + state[Z_VEL]*state[Z_VEL]) 
                      << "\n";
        }
    }
    
    std::cout << "Impact at t=" << t 
              << " x=" << state[X_POS] 
              << " y=" << state[Y_POS] 
              << " z=" << state[Z_POS] << "\n";

    return 0;
}
