/**
 * @file demo_rk45_adaptive.cpp
 * @brief Adaptive RK45 step-size control demonstration.
 *
 * Integrates a high-speed ping-pong ball (k/m = PING_PONG_BALL) under v^2 drag.
 * Drag acceleration scales with v^2, so the launch is the hardest part of the
 * flight and needs a small dt, while the slow terminal-velocity tail allows a
 * much larger one; the constant-Cd drag is smooth, so the embedded error
 * estimate stays well-behaved.
 *
 * Output has two parts: a per-step trace of dt and the error estimate, and a
 * comparison of the adaptive step count against fixed-step RK45 (same order) at
 * matched accuracy.
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/derivative_functions.hpp"
#include "../include/integrators.hpp"
#include "../include/adaptive_driver.hpp"

using namespace StateIndex4D;

// --- Scenario parameters ---
namespace Scenario {
    constexpr double K_OVER_M  = DragRatios::PING_PONG_BALL;  // very high drag (~0.134)
    constexpr double V0        = 200.0;  // m/s launch speed (high, so v^2 drag bites hard)
    constexpr double ANGLE_DEG = 45.0;   // launch angle
    constexpr double H0        = 0.0;    // launch height (m)
    constexpr double DURATION  = 3.0;    // s integration window
    constexpr double DT_INIT   = 0.01;   // initial step-size guess (s)
    constexpr double DT_MAX    = 0.50;   // cap so the calm tail doesn't overstep the window
}

static double speed(const State4D& s) {
    return std::sqrt(s[X_VEL] * s[X_VEL] + s[Y_VEL] * s[Y_VEL]);
}

static State4D launch_state() {
    const double a = Scenario::ANGLE_DEG * M_PI / 180.0;
    return State4D{0.0, Scenario::H0,
                   Scenario::V0 * std::cos(a),
                   Scenario::V0 * std::sin(a)};
}

static double pos_error(const State4D& a, const State4D& b) {
    const double dx = a[X_POS] - b[X_POS];
    const double dy = a[Y_POS] - b[Y_POS];
    return std::sqrt(dx * dx + dy * dy);
}

// Fixed-step RK45 to t=DURATION; reports the number of steps taken.
template <class Deriv>
static State4D run_fixed_to_T(const Deriv& deriv, double dt, long long& steps) {
    State4D s = launch_state();
    double t = 0.0;
    steps = 0;
    while (t < Scenario::DURATION - 1e-12) {
        const double h = std::min(dt, Scenario::DURATION - t);
        s = rk45_step(s, t, h, deriv);
        t += h;
        ++steps;
    }
    return s;
}

// High-order reference solution at t=DURATION (treated as ground truth).
template <class Deriv>
static State4D run_truth(const Deriv& deriv) {
    State4D s = launch_state();
    double t = 0.0;
    const double dt = 1e-4;
    while (t < Scenario::DURATION - 1e-12) {
        const double h = std::min(dt, Scenario::DURATION - t);
        s = rk8_step(s, t, h, deriv);
        t += h;
    }
    return s;
}

struct AdaptiveStats {
    State4D   final_state{};
    long long accepted = 0;
    long long rejected = 0;
    double    min_dt   = std::numeric_limits<double>::infinity();
    double    max_dt   = 0.0;
    double    t_min_dt = 0.0;  // when the smallest step was taken
    double    v_min_dt = 0.0;  // speed there (to tie small dt back to the physics)
};

// Adaptive RK45 integration to t=DURATION. The observer prints the per-step
// trace and records the smallest controller-chosen step (its start-of-step
// time and speed).
template <class Deriv>
static AdaptiveStats run_adaptive(const Deriv& deriv, double tol, bool trace) {
    AdaptiveStats st;
    State4D y = launch_state();

    AdaptiveControl ctrl;
    ctrl.dt_init = Scenario::DT_INIT;
    ctrl.dt_max  = Scenario::DT_MAX;
    ctrl.dt_min  = 1e-9;
    // safety 0.9, max_grow 2.0, min_shrink 0.1, max_shrink_on_reject 0.5 and
    // error_order 4 (exponent 1/5) come from the AdaptiveControl defaults.

    if (trace) {
        std::cout << " step   t(s)    dt(s)      v(m/s)    err est    result\n";
        std::cout << " --------------------------------------------------------------------\n";
    }

    auto observe = [&](const StepInfo<State4D>& info) {
        const double v = speed(info.y);
        if (trace) {
            double err = 0.0;
            for (double e : info.error) err = std::max(err, std::abs(e));
            std::cout << std::setw(5) << info.attempt << "  "
                      << std::fixed << std::setprecision(3) << std::setw(6) << info.t << "  "
                      << std::scientific << std::setprecision(2) << std::setw(8) << info.dt << "  "
                      << std::fixed << std::setprecision(2) << std::setw(8) << v << "   "
                      << std::scientific << std::setprecision(2) << std::setw(8) << err << "   "
                      << (info.accepted ? "accept" : "REJECT")
                      << "\n";
        }
        // Smallest controller-chosen step: exclude the warm-up attempt (#1 uses
        // the raw DT_INIT guess) and the clamped final step; record start-of-step
        // (t, v), where the difficulty that forced this dt lives.
        if (info.accepted && !info.clamped && info.attempt > 1) {
            if (info.dt < st.min_dt) { st.min_dt = info.dt; st.t_min_dt = info.t; st.v_min_dt = v; }
            if (info.dt > st.max_dt) st.max_dt = info.dt;
        }
    };

    const AdaptiveResult res = integrate_adaptive(
        y, 0.0, Scenario::DURATION, ctrl, deriv,
        MaxAbsNorm{tol}, /*event=*/nullptr, observe);

    st.final_state = y;
    st.accepted    = res.accepted;
    st.rejected    = res.rejected;
    return st;
}

int main() {
    ProjectileContext ctx{GRAVITY_EARTH, Scenario::K_OVER_M};
    auto deriv = [&ctx](const State4D& s, double t, State4D& out) {
        drag_deriv_v_squared(s, t, out, ctx);
    };

    std::cout << "=== Adaptive RK45 Demo: High-Speed Ping-Pong Ball (smooth v^2 drag) ===\n\n";
    std::cout << "Scenario: ping-pong ball (k/m = " << std::fixed << std::setprecision(3)
              << Scenario::K_OVER_M << "), v0 = " << std::setprecision(0) << Scenario::V0
              << " m/s, launch " << Scenario::ANGLE_DEG << " deg, smooth v^2 drag.\n"
              << "Integrated to t = " << Scenario::DURATION << " s. Drag accel ~ v^2, so the launch\n"
              << "is violent (small dt) and the slow terminal-velocity tail is gentle (large dt).\n\n";

    // --- 1. Per-step trace ---------------------------------------------------
    const double trace_tol = 1e-4;
    std::cout << "--- Per-step trace (tol = " << std::scientific << std::setprecision(0)
              << trace_tol << ") ---\n";
    const AdaptiveStats tr = run_adaptive(deriv, trace_tol, /*trace=*/true);
    std::cout << " --------------------------------------------------------------------\n";
    std::cout << "Accepted steps: " << tr.accepted << "    Rejected steps: " << tr.rejected << "\n";
    std::cout << "Smallest step:  " << std::scientific << std::setprecision(2) << tr.min_dt
              << " s  (at t = " << std::fixed << std::setprecision(2) << tr.t_min_dt
              << " s, v = " << tr.v_min_dt << " m/s)\n";
    std::cout << "Largest step:   " << std::scientific << std::setprecision(2) << tr.max_dt
              << " s\n\n";

    // --- 2. Adaptive vs. fixed-step efficiency -------------------------------
    const State4D truth = run_truth(deriv);
    std::cout << "--- Why adaptive? Steps needed to reach a given accuracy at t = "
              << std::fixed << std::setprecision(0) << Scenario::DURATION << " s ---\n";
    std::cout << "Reference (RK8, dt = 1e-4): x = " << std::fixed << std::setprecision(4)
              << truth[X_POS] << " m, y = " << truth[Y_POS] << " m\n\n";
    std::cout << "    tol      adaptive steps (acc+rej)   adaptive err (m)   fixed RK45 steps for same err\n";
    std::cout << "  --------------------------------------------------------------------------------------\n";

    for (double tol : {1e-3, 1e-5, 1e-7}) {
        const AdaptiveStats a = run_adaptive(deriv, tol, /*trace=*/false);
        const double aerr = pos_error(a.final_state, truth);
        const long long a_total = a.accepted + a.rejected;

        // Smallest number of *uniform* RK45 steps that matches the adaptive
        // accuracy. Bracket the threshold by doubling the step count, then
        // bisect on the integer count. Searching N directly (dt = DURATION/N)
        // rather than only halving dt avoids quantizing the answer to powers of
        // two. Assumes the error decreases with N, which holds to leading order
        // for a convergent method.
        auto fixed_err = [&](long long N) {
            long long steps = 0;
            const State4D fs = run_fixed_to_T(deriv, Scenario::DURATION / double(N), steps);
            return pos_error(fs, truth);
        };
        long long fixed_steps = -1;
        long long hi = 1;
        // Double until the error is definitively acceptable. Phrased as
        // !(err <= aerr) so that a diverged coarse step (err == NaN/inf, where
        // every comparison is false) counts as "not yet good enough" and we
        // keep increasing N rather than stopping.
        while (hi <= (1LL << 28) && !(fixed_err(hi) <= aerr)) hi <<= 1;
        if (fixed_err(hi) <= aerr) {
            long long lo = hi >> 1;   // largest count still insufficient (0 if hi == 1)
            while (lo + 1 < hi) {
                const long long mid = lo + (hi - lo) / 2;
                if (fixed_err(mid) <= aerr) hi = mid; else lo = mid;
            }
            fixed_steps = hi;
        }

        std::ostringstream acc_rej;
        acc_rej << a_total << " (" << a.accepted << "+" << a.rejected << ")";
        std::cout << "  " << std::scientific << std::setprecision(0) << std::setw(7) << tol
                  << "      " << std::left << std::setw(19) << acc_rej.str() << std::right
                  << std::scientific << std::setprecision(2) << std::setw(10) << aerr
                  << "       " << std::setw(10);
        if (fixed_steps >= 0) std::cout << fixed_steps;
        else                  std::cout << "(> 2^28)";
        std::cout << "\n";
    }
    std::cout << "\nThe fixed-step count is dominated by the violent launch: a uniform dt must be\n"
                 "small enough for the first milliseconds at v0, then wastes that resolution over\n"
                 "the entire calm tail. Adaptive spends its small steps only at the start and\n"
                 "strides through the terminal-velocity descent -- here a ~10x step saving that\n"
                 "widens as the accuracy target tightens. (The adaptive count includes rejected\n"
                 "attempts, which also cost evaluations.)\n";

    return 0;
}
