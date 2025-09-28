#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <functional>
#include <string>
#include "constants.hpp"

//#define DEBUG // define before local includes to enable debug mode in those files
#include "golden_section.hpp"
#include "integrators.hpp"

enum StateIndex { X_POS = 0, Y_POS = 1, X_VEL = 2, Y_VEL = 3 };

void drag_deriv(const State& s, double t, State& deriv, double g, double k_over_m) {
    const double v = sqrt(s[X_VEL]*s[X_VEL] + s[Y_VEL]*s[Y_VEL]);
    deriv[X_POS] = s[X_VEL];  // dx/dt = vx
    deriv[Y_POS] = s[Y_VEL];  // dy/dt = vy
    deriv[X_VEL] = -k_over_m * v * s[X_VEL];  // dvx/dt = -k/m * v * vx
    deriv[Y_VEL] = -g - k_over_m * v * s[Y_VEL];  // dvy/dt = -g -k/m * v * vy
}

void no_drag_deriv(const State& s, double t, State& deriv, double g) {
    deriv[X_POS] = s[X_VEL];  // dx/dt = vx
    deriv[Y_POS] = s[Y_VEL];  // dy/dt = vy
    deriv[X_VEL] = 0.0;       // dvx/dt = 0
    deriv[Y_VEL] = -g;        // dvy/dt = -g
}

struct ScenarioResult {
    double angle;
    double velocity;
    double height;
    double deltaT;
    double distance;
};

ScenarioResult simulate_trajectory(double angle_deg, double v0, double h0, 
                                   double dt, double g, double k_over_m = 0,
                                   IntegratorFunc<State> integrator = euler_step) {                              
    const double angle_rad = angle_deg * M_PI / 180.0;
    State state = {0.0, h0, v0*cos(angle_rad), v0*sin(angle_rad)};
    double t = 0.0;
    double last_x = 0.0, last_y = h0;
    
    auto deriv = [&](const State& s, double t, State& d) {
            drag_deriv(s, t, d, g, k_over_m);
    };

    while (state[1] >= 0.0) {
        last_x = state[0];
        last_y = state[1];
        state = integrator(state, t, dt, deriv);
        t += dt;
    }

    // Linear interpolation for impact point
    const double impact_x = last_x + (0 - last_y) * (state[0] - last_x)/(state[1] - last_y);
    
    return {angle_deg, v0, h0,dt, impact_x};
}



int main() {
    const double g = 9.81;
    const double deltaT = 0.001;
    const double v0 = 50;
    const double h0 = 0.0;
    const double angle_tolerance = 0.000001;

    //const double k_over_m = 0.0;   // No Drag scenario (validation testing)
    // Vacuum 0.0, GolfBall .0025, PingPongBall .01
    const double k_over_m = 0.0025;

    auto distance_func = [&](double angle_deg) {
        auto result = simulate_trajectory(angle_deg, v0, h0, deltaT, g, k_over_m, rk4_step);
        return result.distance;
    };

    // Create a wide initial bracket that definitely contains the maximum
    double a = 10.0;
    double b = 46.0;
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\nRunning golden-section maximization:\n";
    double optimal_angle = golden_section_search_max(distance_func, a, b, angle_tolerance);

    // Final verification with tolerance check
    std::cout << "\nFinal verification near optimum (tolerance = " << angle_tolerance << "°):\n";
    for (double angle = optimal_angle - 1.0; angle <= optimal_angle + 1.0; angle += 0.1) {
        double dist = distance_func(angle);
        double diff = std::abs(angle - optimal_angle);
        
        std::cout << angle << "°: " << dist << " m";
        
        if (diff <= angle_tolerance) {
            std::cout << " <-- WITHIN TOLERANCE OF OPTIMUM";
        }
        std::cout << "\n";
    }

    return 0;
}