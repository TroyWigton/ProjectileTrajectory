/**
 * @file adaptive_driver_test.cpp
 * @brief Tests for the adaptive integration driver (include/adaptive_driver.hpp).
 *
 * Returns non-zero on any failure. Covers accuracy against analytic and RK8
 * references, agreement with a direct inline controller, event localization,
 * and integration over DynamicState (a small N-body system).
 */
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <limits>

#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/integrators.hpp"
#include "../include/adaptive_driver.hpp"
#include "../include/derivative_functions.hpp"

static bool g_all_passed = true;

static void check(bool cond, const std::string& name, const std::string& detail = "") {
    std::cout << (cond ? "  PASSED: " : "  FAILED: ") << name;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << "\n";
    if (!cond) g_all_passed = false;
}

// ---- helpers ----------------------------------------------------------------
static State4D projectile_launch(double v0, double angle_deg, double h0) {
    using namespace StateIndex4D;
    const double a = angle_deg * M_PI / 180.0;
    return State4D{0.0, h0, v0 * std::cos(a), v0 * std::sin(a)};
}

// Minimal 2D N-body gravity over DynamicState; per body layout [x, y, vx, vy].
struct NBodyCtx { int n; std::vector<double> mass; double G; double soft; };

static void nbody_deriv(const DynamicState& s, double, DynamicState& d, const NBodyCtx& c) {
    for (int i = 0; i < c.n; ++i) {
        d[4*i+0] = s[4*i+2];
        d[4*i+1] = s[4*i+3];
        double ax = 0.0, ay = 0.0;
        for (int j = 0; j < c.n; ++j) {
            if (j == i) continue;
            const double dx = s[4*j+0] - s[4*i+0];
            const double dy = s[4*j+1] - s[4*i+1];
            const double r2 = dx*dx + dy*dy + c.soft*c.soft;
            const double inv = 1.0 / (r2 * std::sqrt(r2));
            ax += c.G * c.mass[j] * dx * inv;
            ay += c.G * c.mass[j] * dy * inv;
        }
        d[4*i+2] = ax;
        d[4*i+3] = ay;
    }
}

static double nbody_energy(const DynamicState& s, const NBodyCtx& c) {
    double ke = 0.0, pe = 0.0;
    for (int i = 0; i < c.n; ++i) {
        const double vx = s[4*i+2], vy = s[4*i+3];
        ke += 0.5 * c.mass[i] * (vx*vx + vy*vy);
    }
    for (int i = 0; i < c.n; ++i)
        for (int j = i+1; j < c.n; ++j) {
            const double dx = s[4*j+0] - s[4*i+0];
            const double dy = s[4*j+1] - s[4*i+1];
            const double r = std::sqrt(dx*dx + dy*dy + c.soft*c.soft);
            pe -= c.G * c.mass[i] * c.mass[j] / r;
        }
    return ke + pe;
}

int main() {
    std::cout << std::setprecision(6);

    // --- 1. Accuracy on a smooth scalar ODE: y' = y, y(0)=1 -> e^t ----------
    std::cout << "Test: scaled-norm accuracy on y' = y (truth e^t)\n";
    {
        using S = FixedState<1>;
        auto deriv = [](const S& y, double, S& out) { out[0] = y[0]; };
        S y = {1.0};
        AdaptiveControl ctrl; ctrl.dt_init = 0.01;
        const AdaptiveResult res =
            integrate_adaptive(y, 0.0, 5.0, ctrl, deriv, ScaledNorm{1e-9, 1e-9});
        const double truth = std::exp(5.0);
        const double rel = std::abs(y[0] - truth) / truth;
        check(rel < 1e-6, "final value within 1e-6 of e^5",
              "rel err = " + std::to_string(rel));
        check(res.accepted > 0 && res.t_final > 5.0 - 1e-9, "reached t_end with steps taken");
    }

    // --- 2. Accuracy on real physics: v^2 drag projectile vs RK8 reference ---
    std::cout << "Test: scaled-norm accuracy on v^2 drag projectile (vs RK8)\n";
    {
        ProjectileContext ctx{9.81, 0.0057};
        auto deriv = [&ctx](const State4D& s, double t, State4D& o) {
            drag_deriv_v_squared(s, t, o, ctx);
        };
        const double T = 5.0;
        const State4D y0 = projectile_launch(100.0, 45.0, 0.0);

        // RK8 fine-step reference.
        State4D ref = y0; double t = 0.0; const double dt = 1e-4;
        while (t < T - 1e-12) { const double h = std::min(dt, T - t); ref = rk8_step(ref, t, h, deriv); t += h; }

        State4D y = y0;
        AdaptiveControl ctrl; ctrl.dt_init = 0.01;
        integrate_adaptive(y, 0.0, T, ctrl, deriv, ScaledNorm{1e-10, 1e-10});
        using namespace StateIndex4D;
        const double ex = std::abs(y[X_POS] - ref[X_POS]);
        const double ey = std::abs(y[Y_POS] - ref[Y_POS]);
        check(ex < 1e-4 && ey < 1e-4, "final position within 1e-4 m of RK8 truth",
              "dx=" + std::to_string(ex) + " dy=" + std::to_string(ey));
    }

    // --- 3. Driver (max-abs norm) agrees with a direct inline controller -----
    std::cout << "Test: driver (MaxAbsNorm) matches a direct inline controller\n";
    {
        ProjectileContext ctx{9.81, DragRatios::PING_PONG_BALL};
        auto deriv = [&ctx](const State4D& s, double t, State4D& o) {
            drag_deriv_v_squared(s, t, o, ctx);
        };
        const double T = 3.0, tol = 1e-5, DT_MAX = 0.5;
        const State4D y0 = projectile_launch(200.0, 45.0, 0.0);

        // Independent reference: the same controller written inline.
        State4D yr = y0; double t = 0.0, dt = 0.01; long long acc = 0, rej = 0;
        while (t < T - 1e-12) {
            const double h = std::min(dt, T - t);
            const AdaptiveStepResult<State4D> r = rk45_adaptive_step(yr, t, h, deriv);
            double err = 0.0; for (double e : r.error) err = std::max(err, std::abs(e));
            const double E = err / tol;
            if (E <= 1.0) {
                yr = r.state; t += h; ++acc;
                const double sc = std::min(std::max(0.9*std::pow(1.0/(E+1e-300),0.2),0.1),2.0);
                dt = std::min(h*sc, DT_MAX);
            } else {
                ++rej;
                const double sc = std::min(std::max(0.9*std::pow(1.0/(E+1e-300),0.2),0.1),0.5);
                dt = h*sc;
            }
        }

        State4D y = y0;
        AdaptiveControl ctrl; ctrl.dt_init = 0.01; ctrl.dt_max = DT_MAX; ctrl.dt_min = 1e-9;
        const AdaptiveResult res = integrate_adaptive(y, 0.0, T, ctrl, deriv, MaxAbsNorm{tol});
        using namespace StateIndex4D;
        const double dpos = std::abs(y[X_POS]-yr[X_POS]) + std::abs(y[Y_POS]-yr[Y_POS]);
        check(res.accepted == acc && res.rejected == rej,
              "accepted/rejected counts match inline loop",
              "driver " + std::to_string(res.accepted) + "/" + std::to_string(res.rejected) +
              " vs inline " + std::to_string(acc) + "/" + std::to_string(rej));
        check(dpos < 1e-12, "final state matches inline loop");
    }

    // --- 4. Event localization (MATLAB-style zero crossing) ------------------
    std::cout << "Test: event localization\n";
    {
        // 4a. y' = -1, y0 = 1, event g = y -> crosses 0 at t = 1.
        using S = FixedState<1>;
        auto deriv = [](const S&, double, S& out) { out[0] = -1.0; };
        S y = {1.0};
        EventSpec<S> ev; ev.g = [](double, const S& s) { return s[0]; };
        ev.direction = -1; ev.terminal = true; ev.tol = 1e-12;
        AdaptiveControl ctrl; ctrl.dt_init = 0.1;
        const AdaptiveResult res =
            integrate_adaptive(y, 0.0, 10.0, ctrl, deriv, ScaledNorm{1e-9,1e-9}, &ev);
        check(res.event_terminated, "linear event fired");
        check(std::abs(res.t_final - 1.0) < 1e-8, "event time located at t=1",
              "t_event = " + std::to_string(res.t_final));
        check(std::abs(y[0]) < 1e-8, "state at event ~ 0");
    }
    {
        // 4b. Vacuum projectile ground impact: analytic t = 2 v0 sin(theta)/g.
        ProjectileContext ctx{9.81, 0.0};
        auto deriv = [&ctx](const State4D& s, double t, State4D& o) {
            no_drag_deriv(s, t, o, ctx);
        };
        const double v0 = 100.0, ang = 45.0;
        State4D y = projectile_launch(v0, ang, 0.0);
        using namespace StateIndex4D;
        EventSpec<State4D> ev; ev.g = [](double, const State4D& s) { return s[Y_POS]; };
        ev.direction = -1; ev.terminal = true; ev.tol = 1e-12;
        AdaptiveControl ctrl; ctrl.dt_init = 0.05;
        const AdaptiveResult res =
            integrate_adaptive(y, 0.0, 60.0, ctrl, deriv, ScaledNorm{1e-10,1e-10}, &ev);
        const double t_analytic = 2.0 * v0 * std::sin(ang*M_PI/180.0) / 9.81;
        check(res.event_terminated, "ground-impact event fired");
        check(std::abs(res.t_final - t_analytic) < 1e-4, "impact time matches analytic",
              "t=" + std::to_string(res.t_final) + " vs " + std::to_string(t_analytic));
    }

    // --- 5. State-genericity: 2-body gravity over DynamicState ---------------
    std::cout << "Test: DynamicState (std::vector) N-body integration\n";
    {
        NBodyCtx c; c.n = 2; c.mass = {1.0, 1.0}; c.G = 1.0; c.soft = 1e-3;
        auto deriv = [&c](const DynamicState& s, double t, DynamicState& o) {
            nbody_deriv(s, t, o, c);
        };
        // bodies at (-1,0) and (1,0), opposite transverse velocities (COM fixed)
        DynamicState y = {-1.0, 0.0, 0.0, 0.5,   1.0, 0.0, 0.0, -0.5};
        const double E0 = nbody_energy(y, c);
        AdaptiveControl ctrl; ctrl.dt_init = 0.01;
        const AdaptiveResult res =
            integrate_adaptive(y, 0.0, 2.0, ctrl, deriv, ScaledNorm{1e-10, 1e-10});
        bool finite = true; for (double v : y) if (!std::isfinite(v)) finite = false;
        const double E1 = nbody_energy(y, c);
        const double dE = std::abs((E1 - E0) / E0);
        check(finite && res.accepted > 0, "ran with std::vector state (no NaN, steps taken)");
        check(res.t_final > 2.0 - 1e-9, "reached t_end");
        check(dE < 1e-6, "energy conserved", "rel dE = " + std::to_string(dE));
    }

    std::cout << "\n" << (g_all_passed ? "All adaptive-driver tests passed."
                                       : "SOME TESTS FAILED.") << "\n";
    return g_all_passed ? 0 : 1;
}
