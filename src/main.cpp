#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include "math.h"
#include <functional>
#include <string>

//#define DEBUG

// 4D state: [x_pos, y_pos, x_vel, y_vel]
typedef std::array<double,4> State;

enum StateIndex { X_POS = 0, Y_POS = 1, X_VEL = 2, Y_VEL = 3 };

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

// 8th Order Runge-Kutta
template<typename State, typename DerivativeFunction>
State rk8_step(const State& state, double t, double dt, DerivativeFunction deriv_func) {
    // Dormand-Prince 8(7) Butcher tableau coefficients
    constexpr double a[13][13] = {
        {0},
        {1.0/18.0},
        {1.0/48.0, 1.0/16.0},
        {1.0/32.0, 0, 3.0/32.0},
        {5.0/16.0, 0, -75.0/64.0, 75.0/64.0},
        {3.0/80.0, 0, 0, 3.0/16.0, 3.0/20.0},
        {29443841.0/614563906.0, 0, 0, 77736538.0/692538347.0,
         -28693883.0/1125000000.0, 23124283.0/1800000000.0},
        {16016141.0/946692911.0, 0, 0, 61564180.0/158732637.0,
         22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0},
        {39632708.0/573591083.0, 0, 0, -433636366.0/683701615.0,
         -421739975.0/2616292301.0, 100302831.0/723423059.0,
         790204164.0/839813087.0, 800635310.0/3783071287.0},
        {246121993.0/1340847787.0, 0, 0, -37695042795.0/15268766246.0,
         -309121744.0/1061227803.0, -12992083.0/490766935.0,
         6005943493.0/2108947869.0, 393006217.0/1396673457.0,
         123872331.0/1001029789.0},
        {-1028468189.0/846180014.0, 0, 0, 8478235783.0/508512852.0,
         1311729495.0/1432422823.0, -10304129995.0/1701304382.0,
         -48777925059.0/3047939560.0, 15336726248.0/1032824649.0,
         -45442868181.0/3398467696.0, 3065993473.0/597172653.0},
        {185892177.0/718116043.0, 0, 0, -3185094517.0/667107341.0,
         -477755414.0/1098053517.0, -703635378.0/230739211.0,
         5731566787.0/1027545527.0, 5232866602.0/850066563.0,
         -4093664535.0/808688257.0, 3962137247.0/1805957418.0,
         65686358.0/487910083.0},
        {403863854.0/491063109.0, 0, 0, -5068492393.0/434740067.0,
         -411421997.0/543043805.0, 652783627.0/914296604.0,
         11173962825.0/925320556.0, -13158990841.0/6184727034.0,
         3936647629.0/1978049680.0, -160528059.0/685178525.0,
         248638103.0/1413531060.0}
    };
    
    // Time coefficients c for each stage
    constexpr double c[13] = {
        0.0,
        1.0/18.0,
        1.0/12.0,
        1.0/8.0,
        5.0/16.0,
        3.0/8.0,
        59.0/400.0,
        93.0/200.0,
        5490023248.0/9719169821.0,
        13.0/20.0,
        1201146811.0/1299019798.0,
        1.0,
        1.0
    };
    
    constexpr double b8[13] = {
        14005451.0/335480064.0, 0, 0, 0, 0,
        -59238493.0/1068277825.0, 181606767.0/758867731.0,
        561292985.0/797845732.0, -1041891430.0/1371343529.0,
        760417239.0/1151165299.0, 118820643.0/751138087.0,
        -528747749.0/2220607170.0, 1.0/4.0
    };
    
    State k[13];
    State temp_state;
    
    // Stage 1
    deriv_func(state, t, k[0]);
    
    // Stages 2-13
    for (int stage = 1; stage < 13; ++stage) {
        temp_state = state;
        for (int j = 0; j < stage; ++j) {
            for (size_t i = 0; i < state.size(); ++i) {
                temp_state[i] += dt * a[stage][j] * k[j][i];
            }
        }
        deriv_func(temp_state, t + c[stage] * dt, k[stage]);
    }
    
    // Combine stages for 8th order solution
    State next_state = state;
    for (int stage = 0; stage < 13; ++stage) {
        for (size_t i = 0; i < state.size(); ++i) {
            next_state[i] += dt * b8[stage] * k[stage][i];
        }
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