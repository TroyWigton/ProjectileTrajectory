#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <functional>

//#define DEBUG

// 4D state: [x_pos, y_pos, x_vel, y_vel]
typedef std::array<double,4> State;

template<typename State>
using IntegratorFunc = State(*)(const State&, double, double, 
                               std::function<void(const State&, double, State&)>);

//Forward Euler (Explicit Euler)
template<typename DerivativeFunc>
State euler_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State deriv;
    deriv_func(state, t, deriv);  // Compute derivative at current state
    
    State next_state;
    for (size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + dt * deriv[i];  // Update state
    }
    return next_state;
}

//Heun's Method (Trapezoidal Rule)
template<typename DerivativeFunc>
State heun_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State k1, k2, temp;
    
    // Predictor (Euler step)
    deriv_func(state, t, k1);
    for (size_t i = 0; i < state.size(); ++i) {
        temp[i] = state[i] + dt * k1[i];
    }
    
    // Corrector (Average slopes)
    deriv_func(temp, t + dt, k2);
    State next_state;
    for (size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + 0.5 * dt * (k1[i] + k2[i]);
    }
    return next_state;
}

//Fourth Order Runge-Kutta
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

void drag_deriv(const State& s, double t, State& deriv, double g, double k_over_m) {
    const double v = sqrt(s[2]*s[2] + s[3]*s[3]);
    deriv[0] = s[2];  // dx/dt = vx
    deriv[1] = s[3];  // dy/dt = vy
    deriv[2] = -k_over_m * v * s[2];  // dvx/dt = -k/m * v * vx
    deriv[3] = -g - k_over_m * v * s[3];  // dvy/dt = -g -k/m * v * vy
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
    const double g = 9.81;   // m/s²
    const double k_over_m = 0.0025;  // density * Cd * A / (2m) 
    //const double k_over_m = 0.0;   // No Drag scenario (validation testing)
    // Vacuum 0.0, GolfBall .0025, PingPongBall .01

    const double deltaT = .001;
    const double v0 = 50;
    const double h0 = 0.0;

    std::vector<ScenarioResult> results;
    
    // Parameter sweep
    for (double angle = 40.0; angle <= 45.1; angle += .001) {
        results.push_back(simulate_trajectory(angle, v0, h0, deltaT, g, k_over_m,rk4_step));
    }

    // Find maximum distance scenario
    auto max_drag = *std::max_element(results.begin(), results.end(),
    [](const auto& a, const auto& b) {
        return a.distance < b.distance;
    });

    std::cout << "Using k_over_m = " << k_over_m << std::endl
              <<  "Maximum distance with drag: " 
              <<  max_drag.distance
              << " m at " 
              <<  max_drag.angle 
              << " degrees\n";
    }