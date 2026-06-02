#ifndef ADAPTIVE_DRIVER_HPP
#define ADAPTIVE_DRIVER_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <algorithm>
#include "integrators.hpp"

// Adaptive step-size driver for embedded one-step integrators.
//
// integrate_adaptive() runs a step/accept/reject loop around an embedded
// stepper (default rk45_adaptive_step, which returns {state, error}): it sizes
// each step from the error estimate, rejects and retries when the error exceeds
// tolerance, and can locate zero crossings of an event function g(t, y).
//
// The loop is generic over the state type, and the error measure is a separate
// norm functor, so the same driver handles projectile, N-body, and collision
// problems. The integration method itself lives in integrators.hpp.

namespace adaptive_detail {
    // Clamp x into [lo, hi].
    inline double clamp(double x, double lo, double hi) {
        return std::min(std::max(x, lo), hi);
    }

    // Puts a parameter in a non-deduced context so State is deduced only from
    // `y`; lets callers pass nullptr (no event) or a bare lambda (observer).
    template <class T> struct type_identity { using type = T; };
    template <class T> using nondeduced = typename type_identity<T>::type;
}

// ----- Controller parameters (resize dynamics; tolerance lives in the norm) --
struct AdaptiveControl {
    double dt_init              = 1e-3;   // initial step-size guess
    double dt_min               = 1e-12;  // abort if a rejected step shrinks below this
    double dt_max               = std::numeric_limits<double>::infinity();
    double safety               = 0.9;    // safety factor on the predicted dt
    double max_grow             = 2.0;    // accept: cap dt growth per step
    double min_shrink           = 0.1;    // floor on the resize factor
    double max_shrink_on_reject = 0.5;    // reject: force at least this much shrink
    int    error_order          = 4;      // embedded error order p; exponent = 1/(p+1)
};

// ----- Error-norm policies ---------------------------------------------------
// A norm maps the error estimate to a dimensionless ratio E; a step is accepted
// when E <= 1, and the next dt is sized from E. Any functor with this signature
// works.

// Scaled RMS norm: root-mean-square over components of err_i / (atol + rtol*|y_i|).
// Suits states whose components differ in magnitude (e.g. positions vs velocities).
struct ScaledNorm {
    double atol = 1e-6;
    double rtol = 1e-6;
    template <class State>
    double operator()(const State& err, const State& y, const State& y_new) const {
        const std::size_t n = err.size();
        if (n == 0) return 0.0;
        double acc = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double scale = atol + rtol * std::max(std::abs(y[i]), std::abs(y_new[i]));
            const double r = err[i] / scale;
            acc += r * r;
        }
        return std::sqrt(acc / static_cast<double>(n));
    }
};

// Single absolute tolerance against the infinity-norm of the error: E = max|err|/tol.
struct MaxAbsNorm {
    double tol = 1e-6;
    template <class State>
    double operator()(const State& err, const State& /*y*/, const State& /*y_new*/) const {
        double m = 0.0;
        for (double e : err) m = std::max(m, std::abs(e));
        return m / tol;
    }
};

// ----- Events (MATLAB-style zero-crossing location) --------------------------
template <class State>
struct EventSpec {
    std::function<double(double t, const State& y)> g;  // event function; root = event
    int    direction = 0;      // 0: any crossing, +1: rising only, -1: falling only
    bool   terminal  = true;   // stop integration at the event
    double tol       = 1e-10;  // localization tolerance on the bracketing sub-step width
};

// ----- Per-attempt callback payload ------------------------------------------
template <class State>
struct StepInfo {
    long long    attempt;      // 1-based attempt index (counts rejects too)
    double       t;            // time at the START of the step
    double       dt;           // step size attempted (clamped to t_end if `clamped`)
    const State& y;            // state at the START of the step
    const State& error;        // raw embedded error estimate for this attempt
    double       error_ratio;  // E from the norm (accepted iff <= 1)
    bool         accepted;
    bool         clamped;      // step was clamped to land exactly on t_end
};

// ----- Result / statistics ---------------------------------------------------
struct AdaptiveResult {
    long long accepted = 0;
    long long rejected = 0;
    double    min_dt   = std::numeric_limits<double>::infinity();  // over accepted, unclamped steps
    double    max_dt   = 0.0;
    double    t_final  = 0.0;   // final time reached (event time if terminated early)
    bool      event_terminated = false;
    long long event_index      = -1;  // 0 if the event fired, -1 otherwise
};

namespace adaptive_detail {
    // Default embedded stepper: Dormand-Prince RK45 returning {state, error}.
    struct Rk45Stepper {
        template <class State, class Deriv>
        AdaptiveStepResult<State> operator()(const State& y, double t, double dt, Deriv deriv) const {
            return rk45_adaptive_step(y, t, dt, deriv);
        }
    };

    inline bool is_crossing(double g0, double g1, int direction) {
        if (g0 == 0.0) return false;  // sitting on the event at the step start; don't re-fire
        const bool sign_change = (g0 < 0.0 && g1 >= 0.0) || (g0 > 0.0 && g1 <= 0.0);
        if (!sign_change) return false;
        if (direction > 0) return g1 > g0;   // rising through zero
        if (direction < 0) return g1 < g0;   // falling through zero
        return true;
    }

    // Locate the crossing within (t, t+h] by bisecting the sub-step width,
    // re-stepping from the start state each trial. g0 = g(t, y0). Returns the
    // crossing time and the state there.
    template <class State, class Deriv, class Stepper>
    void localize_event(const State& y0, double t, double h, Deriv deriv, Stepper step,
                        const EventSpec<State>& ev, double g0,
                        double& t_event, State& y_event) {
        double lo = 0.0, hi = h;
        for (int it = 0; it < 200 && (hi - lo) > ev.tol; ++it) {
            const double mid = 0.5 * (lo + hi);
            const State y_mid = step(y0, t, mid, deriv).state;
            const double g_mid = ev.g(t + mid, y_mid);
            const bool mid_opposite = (g0 < 0.0) ? (g_mid >= 0.0) : (g_mid <= 0.0);
            if (mid_opposite) hi = mid; else lo = mid;
        }
        t_event = t + hi;
        y_event = step(y0, t, hi, deriv).state;
    }
}

// ----- The driver ------------------------------------------------------------
// Integrate `y` from `t0` to `t_end` in place, returning statistics. Pass an
// `event` for zero-crossing location and an `observe` callback to inspect each
// attempt.
template <class State, class Deriv,
          class Norm    = ScaledNorm,
          class Stepper = adaptive_detail::Rk45Stepper>
AdaptiveResult integrate_adaptive(
        State& y, double t0, double t_end,
        const AdaptiveControl& ctrl, Deriv deriv,
        Norm norm = Norm{},
        const EventSpec<adaptive_detail::nondeduced<State>>* event = nullptr,
        std::function<void(const StepInfo<adaptive_detail::nondeduced<State>>&)> observe = {},
        Stepper step = Stepper{}) {
    AdaptiveResult res;
    double t  = t0;
    double dt = ctrl.dt_init;
    const double expo = 1.0 / (static_cast<double>(ctrl.error_order) + 1.0);
    long long attempt = 0;

    double g_prev = event ? event->g(t, y) : 0.0;

    while (t < t_end - 1e-12) {
        const double h = std::min(dt, t_end - t);
        const bool clamped = (h < dt);

        const AdaptiveStepResult<State> r = step(y, t, h, deriv);
        const double E = norm(r.error, y, r.state);
        const bool accepted = (E <= 1.0);
        ++attempt;

        if (observe) {
            const StepInfo<State> info{attempt, t, h, y, r.error, E, accepted, clamped};
            observe(info);
        }

        if (accepted) {
            if (event) {
                const double g_new = event->g(t + h, r.state);
                if (adaptive_detail::is_crossing(g_prev, g_new, event->direction)) {
                    double t_event;
                    State  y_event;
                    adaptive_detail::localize_event(y, t, h, deriv, step, *event,
                                                    g_prev, t_event, y_event);
                    y = y_event;
                    t = t_event;
                    res.t_final = t;
                    res.event_terminated = true;
                    res.event_index = 0;
                    if (event->terminal) return res;
                    g_prev = event->g(t, y);  // re-arm at the crossing
                    continue;
                }
                g_prev = g_new;
            }
            y = r.state;
            t += h;
            ++res.accepted;
            if (!clamped) {
                if (h < res.min_dt) res.min_dt = h;
                if (h > res.max_dt) res.max_dt = h;
            }
            const double scale = adaptive_detail::clamp(
                ctrl.safety * std::pow(1.0 / (E + 1e-300), expo),
                ctrl.min_shrink, ctrl.max_grow);
            dt = std::min(h * scale, ctrl.dt_max);
        } else {
            ++res.rejected;
            const double scale = adaptive_detail::clamp(
                ctrl.safety * std::pow(1.0 / (E + 1e-300), expo),
                ctrl.min_shrink, ctrl.max_shrink_on_reject);
            dt = h * scale;
            if (dt < ctrl.dt_min) { res.t_final = t; return res; }  // step-size underflow
        }
    }
    res.t_final = t;
    return res;
}

#endif // ADAPTIVE_DRIVER_HPP
