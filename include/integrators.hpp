#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <array>
#include <functional>

// Represents system dynamics (dy/dt = f(t, y)). Calculates the rate of change for the current state.
template <typename State>
using SystemDerivative = std::function<void(const State& current_state, double current_time, State& derivative_output)>;

// Numerical integration step function. Advances the state by one time step (dt) given the system dynamics.
template <typename State>
using SystemIntegrator = std::function<State(const State& current_state, double current_time, double time_step, SystemDerivative<State> derivative_func)>;

//Forward Euler (Explicit Euler)
template<typename State, typename DerivativeFunc>
State euler_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State deriv;
    deriv_func(state, t, deriv);
    
    State next_state;
    for (size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + dt * deriv[i];
    }
    return next_state;
}


//Heun's Method (Trapezoidal Rule)
template<typename State, typename DerivativeFunc>
State heun_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State k1, k2, temp;
    
    deriv_func(state, t, k1);
    for (size_t i = 0; i < state.size(); ++i) {
        temp[i] = state[i] + dt * k1[i];
    }
    
    deriv_func(temp, t + dt, k2);
    State next_state;
    for (size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + 0.5 * dt * (k1[i] + k2[i]);
    }
    return next_state;
}


//Fourth Order Runge-Kutta
template<typename State, typename DerivativeFunc>
State rk4_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State k1, k2, k3, k4, temp;
    
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

    State next_state;
    for (size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + (dt/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
    return next_state;
}

// RK45 Implementation Details (Shared Coefficients and Logic)
namespace rk45_detail {
    constexpr double c[7] = {0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};
    constexpr double a[7][6] = {
        {0},
        {1.0/5.0},
        {3.0/40.0, 9.0/40.0},
        {44.0/45.0, -56.0/15.0, 32.0/9.0},
        {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0},
        {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0},
        {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0}
    };
    constexpr double b5[7] = {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0};
    constexpr double error_coeffs[7] = {
        35.0/384.0 - 5179.0/57600.0,
        0 - 0,
        500.0/1113.0 - 7571.0/16695.0,
        125.0/192.0 - 393.0/640.0,
        -2187.0/6784.0 - (-92097.0/339200.0),
        11.0/84.0 - 187.0/2100.0,
        0 - 1.0/40.0
    };

    // Computes the 7 intermediate derivatives (stages or 'k' values) for the RK45 step.
    //
    // Algorithm:
    // This evaluates the system derivative function 7 times within the interval [t, t+dt].
    // Each stage uses a weighted sum of previous stages (Explicit RK) to probe the slope
    // at specific points defined by the Butcher tableau. These 7 'k' values capture the
    // curvature of the function and are shared by both the 4th and 5th order solutions.
    template<typename State, typename DerivativeFunc>
    void compute_stages(const State& state, double t, double dt, DerivativeFunc deriv_func, State (&k)[7]) {
        State temp;
        deriv_func(state, t, k[0]);

        for (int i = 1; i < 7; ++i) {
            temp = state;
            for (int j = 0; j < i; ++j) {
                // Determine if a[i][j] is effectively non-zero to avoid unnecessary FLOPs
                if (a[i][j] != 0.0) {
                    for (size_t s = 0; s < state.size(); ++s) {
                        temp[s] += dt * a[i][j] * k[j][s];
                    }
                }
            }
            deriv_func(temp, t + c[i] * dt, k[i]);
        }
    }
}

// Dormand-Prince 5(4) (RK45) - Fixed Step Mode
// Returns the 5th order solution. Can be used as a high-accuracy fixed-step integrator.
//
// How it works:
// Uses the 7 pre-computed stages to construct a 5th-order approximation of the next state.
// Unlike the adaptive version, this function ignores the error estimate coefficients,
// acting purely as a solver that is more accurate per-step than RK4.
template<typename State, typename DerivativeFunc>
State rk45_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State k[7];
    rk45_detail::compute_stages(state, t, dt, deriv_func, k);

    State next_state = state;
    for (size_t s = 0; s < state.size(); ++s) {
        for (int i = 0; i < 7; ++i) {
            if (rk45_detail::b5[i] != 0.0) {
                next_state[s] += dt * rk45_detail::b5[i] * k[i][s];
            }
        }
    }
    return next_state;
}

// Result structure for adaptive integration
template<typename State>
struct AdaptiveStepResult {
    State state; // The 5th order solution
    State error; // The estimated error (difference between 5th and 4th order)
};

// Dormand-Prince 5(4) with Error Estimation - Adaptive Step Mode
// Returns both the high-order (5th) solution and an error estimate.
//
// Role in Adaptive Control:
// The returned 'error' is the difference between the 4th and 5th order solutions.
// This acts as the error signal for a step-size feedback loop:
// - Large Error Estimate -> The step size 'dt' should be reduced to maintain accuracy.
// - Small Error Estimate -> The step size 'dt' can be increased to improve performance.
template<typename State, typename DerivativeFunc>
AdaptiveStepResult<State> rk45_adaptive_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    State k[7];
    rk45_detail::compute_stages(state, t, dt, deriv_func, k);

    AdaptiveStepResult<State> result;
    result.state = state;
    
    // Initialize error with zeros
    for(size_t s=0; s<state.size(); ++s) {
        result.error[s] = 0.0;
    }

    // Combine results to calculate state update and error estimate
    for (size_t s = 0; s < state.size(); ++s) {
        for (int i = 0; i < 7; ++i) {
            if (rk45_detail::b5[i] != 0.0) {
                result.state[s] += dt * rk45_detail::b5[i] * k[i][s];
            }
            if (rk45_detail::error_coeffs[i] != 0.0) {
                result.error[s] += dt * rk45_detail::error_coeffs[i] * k[i][s];
            }
        }
    }
    
    return result;
}

// 8th Order Runge-Kutta
template<typename State, typename DerivativeFunc>
State rk8_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
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

#endif // INTEGRATORS_HPP
