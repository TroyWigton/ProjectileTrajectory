#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <functional>
#include <string>
#include "integrators.hpp"

//#define DEBUG

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

// Modified golden-section search for maximization
template<typename Func>
double golden_section_search_max(Func f, double a, double b, double tol, bool verbose = false) {
    const double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;
    const double inv_golden_ratio = 1.0 / golden_ratio;
    
    double c = b - (b - a) * inv_golden_ratio;
    double d = a + (b - a) * inv_golden_ratio;
    double fc = f(c);
    double fd = f(d);
    
    if (verbose) {
        std::cout << "Starting golden-section maximization in [" << a << ", " << b << "]\n";
        std::cout << "Initial points: " << c << "° (" << fc << "m), " 
                  << d << "° (" << fd << "m)\n";
    }
    
    while (std::abs(c - d) > tol) {
        if (fc > fd) {  // We want to keep the higher value for maximization
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) * inv_golden_ratio;
            fc = f(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) * inv_golden_ratio;
            fd = f(d);
        }
        
        if (verbose) {
            std::cout << "New bracket: [" << a << ", " << b << "], "
                      << "points: " << c << "° (" << fc << "m), "
                      << d << "° (" << fd << "m)\n";
        }
    }
    
    double optimal = (a + b) / 2.0;
    if (verbose) {
        std::cout << "Converged to " << optimal << " degrees with distance " 
                  << f(optimal) << " m\n";
    }
    return optimal;
}

int main() {
    const double g = 9.81;
    const double deltaT = 0.001;
    const double v0 = 50;
    const double h0 = 0.0;
    const double angle_tolerance = 0.0000000001;

    const double k_over_m = 0.0;   // No Drag scenario (validation testing)
    // Vacuum 0.0, GolfBall .0025, PingPongBall .01
    //const double k_over_m = 0.0025;

    auto distance_func = [&](double angle_deg) {
        auto result = simulate_trajectory(angle_deg, v0, h0, deltaT, g, k_over_m, rk4_step);
        return result.distance;
    };

    // Create a wider initial bracket that definitely contains the maximum
    double a = 10.0;
    double b = 46.0;
    
    // First perform a coarse scan to verify unimodality
    std::cout << "Coarse scan to verify maximum is in bracket:\n";
    for (double angle = a; angle <= b; angle += 2.0) {
        double dist = distance_func(angle);
        std::cout << angle << "°: " << dist << " m\n";
    }

    // Run optimization with careful initialization
    std::cout << "\nRunning golden-section maximization:\n";
    double optimal_angle = golden_section_search_max(distance_func, a, b, angle_tolerance, true);

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