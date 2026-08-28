#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <array>
#include <cstddef>
#include <functional>

// Represents system dynamics (dy/dt = f(t, y)). Calculates the rate of change for the current state.
template <typename State>
using SystemDerivative = std::function<void(const State& current_state, double current_time, State& derivative_output)>;

// Numerical integration step function. Advances the state by one time step (dt) given the system dynamics.
template <typename State>
using SystemIntegrator = std::function<State(const State& current_state, double current_time, double time_step, SystemDerivative<State> derivative_func)>;

// Stage temporaries (k1, k2, ..., temp) are declared thread_local inside each
// integrator: one set per thread per template instantiation, lazily sized from
// the input state on the first call, then reused across subsequent calls.
//
// For runtime-sized containers (std::vector<double>) this avoids per-step heap
// allocation of the stage workspace — substantial savings when many steps run.
// For compile-time-sized containers (std::array<double, N>) the size check
// folds to a compile-time constant and the conditional initialization is
// elided.
//
// Return values (next_state, result.state, result.error) are plain locals so
// that RVO/NRVO can construct the caller's destination directly, avoiding an
// extra copy on the return path.
//
// Constraint: thread_local workspace is shared across all invocations on the
// same thread. The integrators must NOT be invoked recursively (e.g., a
// derivative function calling the same integrator) — that would clobber the
// in-flight workspace. Standard numerical integration patterns never do this.

//Forward Euler (Explicit Euler)
template<typename State, typename DerivativeFunc>
State euler_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    thread_local State deriv;
    if (deriv.size() != state.size()) { deriv = state; }

    deriv_func(state, t, deriv);

    State next_state = state;
    for (std::size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + dt * deriv[i];
    }
    return next_state;
}

// Stateless functor wrapping euler_step so a chosen method can be passed as a
// compile-time Integrator template argument (e.g. to Simulation) instead of
// through std::function type erasure. The call operator is itself a template,
// so State/DerivativeFunc are deduced fresh at each call site and the whole
// chain inlines. Same pattern as adaptive_driver.hpp's Rk45Stepper.
struct EulerStepper {
    template<typename State, typename DerivativeFunc>
    State operator()(const State& state, double t, double dt, DerivativeFunc deriv_func) const {
        return euler_step(state, t, dt, deriv_func);
    }
};


//Heun's Method (Trapezoidal Rule)
template<typename State, typename DerivativeFunc>
State heun_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    thread_local State k1, k2, temp;
    if (k1.size() != state.size()) { k1 = state; k2 = state; temp = state; }

    deriv_func(state, t, k1);
    for (std::size_t i = 0; i < state.size(); ++i) {
        temp[i] = state[i] + dt * k1[i];
    }

    deriv_func(temp, t + dt, k2);
    State next_state = state;
    for (std::size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + 0.5 * dt * (k1[i] + k2[i]);
    }
    return next_state;
}

// See EulerStepper above for why this exists.
struct HeunStepper {
    template<typename State, typename DerivativeFunc>
    State operator()(const State& state, double t, double dt, DerivativeFunc deriv_func) const {
        return heun_step(state, t, dt, deriv_func);
    }
};


//Fourth Order Runge-Kutta
// The "Classical" Runge-Kutta method.
// A 4th-order method, meaning the local error per step is O(dt^5) and global error is O(dt^4).
//
// Geometric Interpretation of the 4 Stages:
// k1: Slope at the start of the interval (t). This is the Euler predictor.
// k2: Slope at the midpoint (t + 0.5*dt). The state is estimated using k1. This estimates the slope at the center.
// k3: Slope at the midpoint (t + 0.5*dt). The state is re-estimated using k2. This refines the midpoint slope.
// k4: Slope at the end of the interval (t + dt). The state is estimated using k3.
//
// The final state update uses a weighted average of these slopes: (k1 + 2*k2 + 2*k3 + k4) / 6.
// The middle slopes (k2, k3) are given double weight, similar to Simpson's rule for integration,
// which cancels lower-order error terms.
template<typename State, typename DerivativeFunc>
State rk4_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    thread_local State k1, k2, k3, k4, temp;
    if (k1.size() != state.size()) { k1 = state; k2 = state; k3 = state; k4 = state; temp = state; }

    // Stage 1: Derivative at the start
    deriv_func(state, t, k1);
    for (std::size_t i = 0; i < state.size(); ++i)
        temp[i] = state[i] + 0.5 * dt * k1[i];

    // Stage 2: Derivative at midpoint (based on k1)
    deriv_func(temp, t + 0.5*dt, k2);
    for (std::size_t i = 0; i < state.size(); ++i)
        temp[i] = state[i] + 0.5 * dt * k2[i];

    // Stage 3: Derivative at midpoint (based on k2)
    // Note: We sample at the *same time* (t + 0.5*dt) as Stage 2.
    // The difference is that we use k2 (the improved slope estimate) to reach this point, rather than k1.
    deriv_func(temp, t + 0.5*dt, k3);
    for (std::size_t i = 0; i < state.size(); ++i)
        temp[i] = state[i] + dt * k3[i];

    // Stage 4: Derivative at end (based on k3)
    // We estimate the state at the end of the interval (t + dt) by projecting
    // from the start using the refined midpoint slope (k3) across the full step dt.
    deriv_func(temp, t + dt, k4);

    // Final State Update
    // We compute the true next state by taking a weighted average of the four slopes.
    // k1, k4 (endpoints) get weight 1.
    // k2, k3 (midpoints) get weight 2.
    State next_state = state;
    for (std::size_t i = 0; i < state.size(); ++i) {
        next_state[i] = state[i] + (dt/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
    return next_state;
}

// See EulerStepper above for why this exists.
struct RK4Stepper {
    template<typename State, typename DerivativeFunc>
    State operator()(const State& state, double t, double dt, DerivativeFunc deriv_func) const {
        return rk4_step(state, t, dt, deriv_func);
    }
};

// RK45 (Dormand-Prince) Implementation Details
//
// About the Method:
// Dormand-Prince is a member of the Runge-Kutta family of ODE solvers.
// It is an "Embedded" method, meaning it generates two results per step:
// 1. A 5th-order approximation (high accuracy).
// 2. A 4th-order approximation (lower accuracy).
// The difference between these two approximations provides an estimate of the local truncation error.
//
// These coefficients (The Butcher Tableau) satisfy specific algebraic conditions so that
// the Taylor series expansion of the approximate solution matches the true solution up to
// order 5 (for the main result) and order 4 (for the error estimator).
namespace rk45_detail {
    // Stage Time Nodes (c coefficients):
    // Determines the time nodes t_i = t + c_i * dt where we evaluate the system's differential equation.
    //
    // The 'k' values are stage slopes (derivatives dy/dt): each stage samples the
    // derivative function f(t, y) at the time node t + c_i*dt listed below.
    //
    // k1 (c=0)    : Slope at start of step (t).
    // k2 (c=1/5)  : Slope at 20% through the time step (t + 0.2*dt).
    // k3 (c=3/10) : Slope at 30% through the time step.
    // k4 (c=4/5)  : Slope at 80% through the time step.
    // k5 (c=8/9)  : Slope at ~89% through the time step.
    // k6 (c=1)    : Slope at the end of the time step (100%).
    // k7 (c=1)    : Slope at the end (100%). Enables the FSAL property (b5 equals the
    //               last row of 'a'): an FSAL driver could reuse k7 as k1 of the next
    //               step, but this implementation recomputes all stages on every call.
    constexpr double c[7] = {0, 1.0/5.0, 3.0/10.0, 4.0/5.0, 8.0/9.0, 1.0, 1.0};

    // Runge-Kutta Matrix (a coefficients):
    // These define how previous stages contribute to the predictor for the current stage.
    // For an explicit method, the matrix is lower triangular (a[i][j] = 0 for j >= i).
    //
    // Interpretation of indices a[i][j]:
    // - i (Row): The "Target" stage we are currently calculating.
    // - j (Column): The "Source" stage (a previous slope) we are adding to the mix.
    //
    // The value a[i][j] is the weight applied to slope k[j] when calculating the input state for k[i].
    // Example: To calculate stage k[2] (i=2), we mix in specific amounts of k[0] (j=0) and k[1] (j=1).
    // This allows the algorithm to "look back" at previous curvature estimates to refine the current one.
    constexpr double a[7][6] = {
        {0}, // k1 depends on nothing
        {1.0/5.0}, // k2 depends on k1
        {3.0/40.0, 9.0/40.0}, // k3 depends on k1, k2
        {44.0/45.0, -56.0/15.0, 32.0/9.0}, // ...
        {19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0},
        {9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0},
        {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0}
    };

    // Weights (b5 coefficients) for the 5th-order solution:
    // The final result y_{n+1} is a weighted sum of the stages: y_{n+1} = y_n + dt * Σ(b5_i * k_i).
    // Note that b5[1] is 0, meaning the k2 stage does not contribute directly to the 5th-order solution,
    // though it influences it indirectly through subsequent stages k3..k7.
    // Also, b5 matches the last row of 'a' (the predictor for k7), which is the FSAL property.
    constexpr double b5[7] = {35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0};

    // Error Estimate Coefficients (E = b5 - b4):
    // Instead of explicitly computing the 4th-order solution, we compute the difference directly.
    // Error = y_5th - y_4th = sum((b5_i - b4_i) * k_i).
    // This vector represents that difference.
    // Notice the last element is non-zero: b5[6] = 0 but b4[6] = 1/40, so the 4th-order
    // solution (and hence the error estimate) depends on k7 even though the 5th-order result does not.
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
        thread_local State temp;
        if (temp.size() != state.size()) { temp = state; }

        deriv_func(state, t, k[0]);

        for (int i = 1; i < 7; ++i) {
            temp = state;
            for (int j = 0; j < i; ++j) {
                // Determine if a[i][j] is effectively non-zero to avoid unnecessary FLOPs
                if (a[i][j] != 0.0) {
                    for (std::size_t s = 0; s < state.size(); ++s) {
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
// By discarding the error estimate, this function behaves effectively as a standard
// fixed-step integrator, but with 5th-order accuracy (7 stage evaluations per step).
template<typename State, typename DerivativeFunc>
State rk45_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    thread_local State k[7];
    if (k[0].size() != state.size()) { for (auto& ki : k) ki = state; }

    rk45_detail::compute_stages(state, t, dt, deriv_func, k);

    State next_state = state;
    for (std::size_t s = 0; s < state.size(); ++s) {
        for (int i = 0; i < 7; ++i) {
            // b5[i] are the weights for the 5th order solution
            if (rk45_detail::b5[i] != 0.0) {
                next_state[s] += dt * rk45_detail::b5[i] * k[i][s];
            }
        }
    }
    return next_state;
}

// See EulerStepper above for why this exists. Wraps the fixed-step (non
// error-estimating) form of RK45.
struct RK45Stepper {
    template<typename State, typename DerivativeFunc>
    State operator()(const State& state, double t, double dt, DerivativeFunc deriv_func) const {
        return rk45_step(state, t, dt, deriv_func);
    }
};

// Result structure for adaptive integration
template<typename State>
struct AdaptiveStepResult {
    State state; // The 5th order solution
    State error; // The estimated error (difference between 5th and 4th order)
};

// Dormand-Prince 5(4) with Error Estimation - Adaptive Step Mode
// This is the core component of an adaptive ODE solver (like MATLAB's ode45).
//
// Principle:
// Both solutions are produced from the same set of stage evaluations:
// - y_5th: The proposed new state.
// - y_4th: A slightly less accurate state.
//
// The difference (error = y_5th - y_4th) is a computable number that scales with the
// local step size. If this error is:
// - Too large (> tolerance): The step must be rejected and retried with smaller dt.
// - Too small (<< tolerance): The next step size can be aggressively increased.
//
// This lets the solver take small steps at sharp turns (high curvature) and large
// steps on straightaways (linear dynamics).
template<typename State, typename DerivativeFunc>
AdaptiveStepResult<State> rk45_adaptive_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    thread_local State k[7];
    if (k[0].size() != state.size()) { for (auto& ki : k) ki = state; }

    rk45_detail::compute_stages(state, t, dt, deriv_func, k);

    // Return-value members sized from the input state. `result.state` is then
    // accumulated into by the 5th-order weighting loop below; `result.error`
    // is explicitly zeroed first so the difference-coeffs loop accumulates
    // cleanly.
    AdaptiveStepResult<State> result;
    result.state = state;
    result.error = state;

    for(std::size_t s=0; s<state.size(); ++s) {
        result.error[s] = 0.0;
    }

    // Combine results to calculate state update and error estimate
    for (std::size_t s = 0; s < state.size(); ++s) {
        for (int i = 0; i < 7; ++i) {
            // Accumulate 5th order solution
            if (rk45_detail::b5[i] != 0.0) {
                result.state[s] += dt * rk45_detail::b5[i] * k[i][s];
            }
            // Accumulate error difference
            if (rk45_detail::error_coeffs[i] != 0.0) {
                result.error[s] += dt * rk45_detail::error_coeffs[i] * k[i][s];
            }
        }
    }
    
    return result;
}

// 8th Order Runge-Kutta (Prince-Dormand RK8(7)13M)
// These are the RK8(7)13M coefficients of Prince & Dormand (1981): a 13-stage, 8th-order
// formula with an embedded 7th-order solution for error estimation ("13" = stage count,
// "M" = their minimized-truncation-error variant). This implementation is fixed-step --
// it propagates the 8th-order solution only and omits the embedded 7th-order estimator.
//
// Usage:
// Suited for problems with very stringent error tolerances (e.g., 10^-9 to 10^-13)
// or extremely smooth functions (like celestial mechanics/orbits).
// For looser tolerances, the overhead of 13 stages per step often makes it slower than RK4 or RK45.
//
// Reference for these coefficients: P. J. Prince and J. R. Dormand, "High order embedded
// Runge-Kutta formulae," J. Comput. Appl. Math. 7(1), 67-75 (1981).
template<typename State, typename DerivativeFunc>
State rk8_step(const State& state, double t, double dt, DerivativeFunc deriv_func) {
    // Butcher Tableau Matrix A (13 stages)
    // Determines the linear combinations of previous k values to predict the next stage.
    // The exact fractions are derived to minimize the error terms of the 8th order Taylor series.
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

    // Stage Time Nodes (c coefficients):
    // Time offsets t_i = t + c[i]*dt for the 13 stages of the RK8 method.
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
    
    // Weights (b coefficients) for the 8th-order solution:
    // Combined to form the final result from the 13 stages: y_{n+1} = y_n + dt * sum(b8_i * k_i).
    constexpr double b8[13] = {
        14005451.0/335480064.0, 0, 0, 0, 0,
        -59238493.0/1068277825.0, 181606767.0/758867731.0,
        561292985.0/797845732.0, -1041891430.0/1371343529.0,
        760417239.0/1151165299.0, 118820643.0/751138087.0,
        -528747749.0/2220607170.0, 1.0/4.0
    };
    
    thread_local State k[13];
    thread_local State temp_state;
    if (k[0].size() != state.size()) {
        for (auto& ki : k) ki = state;
        temp_state = state;
    }

    // Stage 1
    deriv_func(state, t, k[0]);
    
    // Stages 2-13
    for (int stage = 1; stage < 13; ++stage) {
        temp_state = state;
        for (int j = 0; j < stage; ++j) {
            for (std::size_t i = 0; i < state.size(); ++i) {
                temp_state[i] += dt * a[stage][j] * k[j][i];
            }
        }
        deriv_func(temp_state, t + c[stage] * dt, k[stage]);
    }
    
    // Combine stages for 8th order solution
    State next_state = state;
    for (int stage = 0; stage < 13; ++stage) {
        for (std::size_t i = 0; i < state.size(); ++i) {
            next_state[i] += dt * b8[stage] * k[stage][i];
        }
    }
    
    return next_state;
}

// See EulerStepper above for why this exists.
struct RK8Stepper {
    template<typename State, typename DerivativeFunc>
    State operator()(const State& state, double t, double dt, DerivativeFunc deriv_func) const {
        return rk8_step(state, t, dt, deriv_func);
    }
};

#endif // INTEGRATORS_HPP
