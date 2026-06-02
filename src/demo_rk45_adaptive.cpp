/**
 * @file demo_rk45_adaptive.cpp
 * @brief Demonstration of adaptive Runge-Kutta 4(5) integration.
 *
 * The point of an adaptive stepper is that it places its work where the work
 * is needed: it SHRINKS dt where the dynamics are hard and GROWS dt where they
 * are smooth, all driven by the embedded 4th/5th-order error estimate.
 *
 * Scenario: a ping-pong ball (very high drag, k/m = PING_PONG_BALL) launched at
 * high speed under smooth v^2 drag. Because drag acceleration scales with v^2,
 * the opening of the flight decelerates hard -- at launch the ball sheds speed
 * at thousands of m/s^2 (~5,400 m/s^2) -- so the solution has high curvature
 * there and demands a small dt. As the ball slows toward its gentle
 * terminal-velocity descent the dynamics relax and dt can grow by over an order
 * of magnitude. The drag model is
 * smooth (constant Cd), so the error estimate is well-behaved -- this isolates
 * adaptivity from any non-smoothness.
 *
 * Two lessons in the output:
 *   1. A per-step trace, so you can watch dt start tiny in the high-deceleration opening
 *      and stretch out as the motion calms.
 *   2. A head-to-head: how many *fixed*-step RK45 steps (same method, same
 *      order) are required to match the adaptive run's accuracy. Adaptive wins
 *      because a fixed dt must be small *everywhere* to survive the opening,
 *      while the difficulty is concentrated at the start.
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

using namespace StateIndex4D;

// --- Scenario (tweak these to explore the controller's behavior) ---
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

// Adaptive RK45 to t=DURATION. With trace=true, logs every attempt so the
// reader can watch dt and the error estimate respond to the dynamics.
template <class Deriv>
static AdaptiveStats run_adaptive(const Deriv& deriv, double tol, bool trace) {
    AdaptiveStats st;
    State4D s = launch_state();
    double t  = 0.0;
    double dt = Scenario::DT_INIT;

    if (trace) {
        std::cout << " step   t(s)    dt(s)      v(m/s)    err est    result\n";
        std::cout << " --------------------------------------------------------------------\n";
    }

    long long attempt = 0;
    while (t < Scenario::DURATION - 1e-12) {
        // Clamp the final step so we land exactly on t=DURATION. A clamped step
        // is not representative of the controller's choice, so it is excluded
        // from the dt statistics below.
        const double natural_dt = dt;
        const double h = std::min(dt, Scenario::DURATION - t);
        const bool   clamped = (h < natural_dt);

        const AdaptiveStepResult<State4D> r = rk45_adaptive_step(s, t, h, deriv);

        double err = 0.0;
        for (double e : r.error) err = std::max(err, std::abs(e));

        const bool   accept = (err <= tol);
        const double v      = speed(s);
        ++attempt;

        if (trace) {
            std::cout << std::setw(5) << attempt << "  "
                      << std::fixed << std::setprecision(3) << std::setw(6) << t << "  "
                      << std::scientific << std::setprecision(2) << std::setw(8) << h << "  "
                      << std::fixed << std::setprecision(2) << std::setw(8) << v << "   "
                      << std::scientific << std::setprecision(2) << std::setw(8) << err << "   "
                      << (accept ? "accept" : "REJECT")
                      << "\n";
        }

        if (accept) {
            // Record dt stats against the step's STARTING (t, v) -- that is
            // where the difficulty that forced this dt lives. Exclude the
            // warm-up attempt (attempt == 1 uses the raw DT_INIT guess rather
            // than a controller decision) and any step clamped to land on
            // t=DURATION. (Once the first attempt is rejected, every subsequent
            // dt is controller-derived, so this keeps the genuine smallest step
            // at the high-deceleration launch.)
            if (!clamped && attempt > 1) {
                if (h < st.min_dt) { st.min_dt = h; st.t_min_dt = t; st.v_min_dt = v; }
                if (h > st.max_dt) st.max_dt = h;
            }
            s = r.state;
            t += h;
            ++st.accepted;
            // Step-size controller: target dt_new = dt * safety * (tol/err)^(1/5).
            const double scale = std::min(std::max(
                0.9 * std::pow(tol / (err + 1e-30), 0.2), 0.1), 2.0);
            dt = std::min(h * scale, Scenario::DT_MAX);
        } else {
            ++st.rejected;
            // Reject: shrink and retry the same t (do not advance). Use the
            // error-predicted factor, but force at least a halving so that an
            // over-large dt (e.g. the initial guess fed into the high-deceleration launch)
            // collapses quickly instead of creeping down over many attempts.
            double scale = 0.9 * std::pow(tol / (err + 1e-30), 0.2);
            scale = std::min(scale, 0.5);   // shrink by at least half
            scale = std::max(scale, 0.1);   // but don't overshrink in one go
            dt = h * scale;
            if (dt < 1e-9) {
                std::cerr << "Error: step-size underflow; aborting.\n";
                break;
            }
        }
    }
    st.final_state = s;
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
              << "decelerates hard (~5,400 m/s^2, small dt) and the slow terminal-velocity tail is\n"
              << "gentle (large dt).\n\n";

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
    std::cout << "\ntol is the LOCAL error bound (max-norm over the state vector) -- the embedded\n"
                 "5th-minus-4th-order difference within a single RK45 step. Each step is accepted\n"
                 "when that estimate is <= tol, else rejected and retried with a smaller dt.\n"
                 "adaptive err is the GLOBAL position error at t=3 s vs the RK8 reference.\n";
    std::cout << "\nThe fixed-step count is dominated by the high deceleration at launch: a uniform\n"
                 "dt must be small enough for the first milliseconds at v0, then wastes that\n"
                 "resolution over the entire calm tail. Adaptive spends its small steps only at\n"
                 "the start and strides through the terminal-velocity descent -- here a ~10x step\n"
                 "saving that widens as the accuracy target tightens. (The adaptive count includes\n"
                 "rejected attempts, which also cost evaluations.)\n";

    return 0;
}
