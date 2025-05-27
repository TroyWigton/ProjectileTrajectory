#include <vector>
#include <array>
#include <iostream>
#include <iomanip>  // For std::fixed and std::setprecision
#include "math.h"

//#define DEBUG

// 4D state: [x_pos, y_pos, x_vel, y_vel]
using State = std::array<double, 4>;

template<typename DerivativeFunc>
State rk4_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State k1, k2, k3, k4, temp;
    
    // RK4 stages
    deriv_func(state, t, k1);
    for (size_t i = 0; i < state.size(); ++i)
        temp[i] = state[i] + 0.5 * dt * k1[i];
    
    deriv_func(temp, t + 0.5*dt, k2);
    for (size_t i = 0; i < state.size(); ++i)
        temp[i] = state[i] + 0.5 * dt * k2[i];
    
    deriv_func(temp, t + 0.5*dt, k3);
    for (size_t i = 0; i < state.size(); ++i)
        temp[i] = state[i] + dt * k3[i];
    
    deriv_func(temp, t + dt, k4);

    // Combine results
    State next_state;
    for (size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + (dt/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
    return next_state;
}

void no_drag_deriv(const State& s, double t, State& deriv, double g) {
    deriv[0] = s[2];  // dx/dt = vx
    deriv[1] = s[3];  // dy/dt = vy
    deriv[2] = 0.0;  // dvx/dt = 0
    deriv[3] = -g;    // dvy/dt = -g
}

void drag_deriv(const State& s, double t, State& deriv, double g, double k_over_m) {
    const double v = sqrt(s[2]*s[2] + s[3]*s[3]);
    deriv[0] = s[2];  // dx/dt = vx
    deriv[1] = s[3];  // dy/dt = vy
    deriv[2] = -k_over_m * v * s[2];  // dvx/dt = -k/m * v * vx
    deriv[3] = -g - k_over_m * v * s[3];  // dvy/dt = -g -k/m * v * vy
    
    #ifdef DEBUG
    std::cout << "Drag calculation active (v=" 
    << sqrt(s[2]*s[2] + s[3]*s[3]) 
    << " m/s)\n";
    #endif
}

struct ScenarioResult {
    double angle;
    double velocity;
    double height;
    bool use_drag;
    double deltaT;
    double distance;
};

ScenarioResult simulate_trajectory(double angle_deg, double v0, double h0, 
                                  bool use_drag, double dt, double g, double k_over_m = 0) {
    const double angle_rad = angle_deg * M_PI / 180.0;
    State state = {0.0, h0, v0*cos(angle_rad), v0*sin(angle_rad)};
    double t = 0.0;
    double last_x = 0.0, last_y = h0;
    
    auto deriv = [&](const State& s, double t, State& d) {
        if (use_drag) {
            drag_deriv(s, t, d, g, k_over_m);
        } else {
            no_drag_deriv(s, t, d, g);
        }
    };

    while (state[1] >= 0.0) {
        last_x = state[0];
        last_y = state[1];
        state = rk4_step(state, t, dt, deriv);
        t += dt;
    }

    // Linear interpolation for impact point
    const double impact_x = last_x + (0 - last_y) * (state[0] - last_x)/(state[1] - last_y);
    
    return {angle_deg, v0, h0, use_drag, dt, impact_x};
}

int main() {
    const double g = 9.81;   // m/s²
    const double k_over_m = 0.01;  // Drag coefficient/mass
    const double deltaT = .001;
    const double v0 = 50;
    const double h0 = 0.0;

    std::vector<ScenarioResult> results;
    
    // Parameter sweep
    for (double angle = 35.0; angle <= 45.1; angle += .001) {
        // Simulate both scenarios
        results.push_back(simulate_trajectory(angle, v0, h0, false, deltaT, g));
        results.push_back(simulate_trajectory(angle, v0, h0, true, deltaT, g, k_over_m));
    }

// Find maximum distance scenarios
// For no-drag maximum
auto max_no_drag = *std::max_element(results.begin(), results.end(),
[](const auto& a, const auto& b) {
    return !a.use_drag && a.distance < b.distance;
});

// For drag maximum (corrected comparator)
auto max_drag = *std::max_element(results.begin(), results.end(),
[](const auto& a, const auto& b) {
    if (a.use_drag && b.use_drag) return a.distance < b.distance;
    return !a.use_drag && b.use_drag;
});

    // Output comparison
std::cout << std::fixed << std::setprecision(4);  // Applies to ALL subsequent floats
std::cout << "Maximum distance without drag: " 
          << std::defaultfloat << max_no_drag.distance  // Reset for distance
          << " m at " 
          << std::fixed << max_no_drag.angle 
          << " degrees\n";

std::cout << "Maximum distance with drag: " 
          << std::defaultfloat << max_drag.distance  // Reset for distance
          << " m at " 
          << std::fixed << max_drag.angle 
          << " degrees\n";
}