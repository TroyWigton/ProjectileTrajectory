/**
 * @file convergence_study.cpp
 * @brief Empirical order-of-convergence measurement for the fixed-step integrators.
 *
 * For each integrator, runs a smooth v^2-drag projectile to a fixed time T at a
 * sequence of progressively halved dt values, and reports the global position
 * error at t=T against an ultra-high-precision RK8 reference. For a method of
 * theoretical order p, halving dt should shrink the global error by 2^p — so
 * the empirical order log2(error_prev / error_curr) should approach p.
 *
 * Once the error reaches the precision floor (~1e-10 for this problem, set by
 * the reference's own accumulated rounding), shrinking dt further stops reducing
 * the error, and the empirical order departs from the theoretical value. That
 * breakdown is itself part of the convergence story for high-order methods on
 * smooth problems and is left visible in the output rather than hidden.
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include "constants.hpp"
#include "types.hpp"
#include "integrators.hpp"
#include "derivative_functions.hpp"

// Integrate from initial state s0 over [0, t_end] using the given fixed-step
// integrator at step size dt. Recomputes actual_dt = t_end / round(t_end/dt)
// to land exactly on t_end for fair error comparison.
template<typename Integrator, typename Deriv>
State4D run_to_t_end(Integrator integrator, double t_end, double dt,
                     const State4D& s0, Deriv deriv) {
    State4D state = s0;
    const long long steps = std::lround(t_end / dt);
    const double actual_dt = t_end / static_cast<double>(steps);
    double t = 0.0;
    for (long long i = 0; i < steps; ++i) {
        state = integrator(state, t, actual_dt, deriv);
        t += actual_dt;
    }
    return state;
}

// Euclidean distance between two 2D position vectors carried in the state.
double position_error(const State4D& a, const State4D& b) {
    using namespace StateIndex4D;
    const double dx = a[X_POS] - b[X_POS];
    const double dy = a[Y_POS] - b[Y_POS];
    return std::sqrt(dx*dx + dy*dy);
}

int main() {
    using namespace StateIndex4D;

    // Test scenario: smooth v^2 drag at a moderate launch velocity. Chosen so
    // every integrator can complete the dt sweep in well under a second and
    // the dynamics are smooth (no Cd kinks, no transonic transitions) so the
    // theoretical orders are recoverable.
    const double v0           = 100.0;
    const double launch_deg   = 45.0;
    const double launch_rad   = launch_deg * M_PI / 180.0;
    const double t_end        = 5.0;
    const ProjectileContext ctx{GRAVITY_EARTH, DragRatios::GOLF_BALL};

    auto deriv = [&ctx](const State4D& s, double t, State4D& d) {
        drag_deriv_v_squared(s, t, d, ctx);
    };

    const State4D s0 = {0.0, 0.0,
                        v0 * std::cos(launch_rad),
                        v0 * std::sin(launch_rad)};

    // Ground truth: RK8 at very small dt. 1e-6 over 5s = 5M steps × 13 stages
    // = ~65M derivative evaluations — runs in a few seconds, well below the
    // floating-point precision floor.
    const double ref_dt = 1.0e-6;
    const State4D truth = run_to_t_end(
        rk8_step<State4D, decltype(deriv)>, t_end, ref_dt, s0, deriv);

    std::cout << "Numerical Integrator Convergence Test\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Test problem (smooth v^2 drag):\n";
    std::cout << "  v0           = " << v0 << " m/s\n";
    std::cout << "  launch angle = " << launch_deg << " deg\n";
    std::cout << "  k/m          = " << ctx.k_over_m << " (golf ball)\n";
    std::cout << "  duration T   = " << t_end << " s\n";
    std::cout << "Reference: RK8 at dt = " << ref_dt
              << " (truth = (" << truth[X_POS] << ", " << truth[Y_POS] << ") m)\n";
    std::cout << "\n";
    std::cout << "Each integrator is run at progressively halved dt values; the\n";
    std::cout << "global position error at t=T is computed against the reference.\n";
    std::cout << "Empirical order = log2(error_prev / error_curr). For an order-p\n";
    std::cout << "method, the empirical order should approach p until the error\n";
    std::cout << "reaches the precision floor (~1e-10 for this problem).\n";
    std::cout << "\n";

    // Each integrator gets a dt sequence chosen to keep error well above the
    // floating-point precision floor (~1e-10 for this problem) while still
    // showing several halvings of meaningful convergence behavior. Higher-order
    // methods need larger dt values to keep error in the observable range;
    // lower-order methods need smaller dt to keep error from being uselessly
    // huge.
    const std::vector<double> dts_low_order = {0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125};
    const std::vector<double> dts_rk4       = {0.1, 0.05, 0.025, 0.0125, 0.00625};
    const std::vector<double> dts_rk45      = {0.5, 0.25, 0.125, 0.0625, 0.03125};
    const std::vector<double> dts_rk8       = {1.0, 0.5, 0.25, 0.125, 0.0625};

    auto run_convergence_table = [&](const std::string& name,
                                     int theoretical_order,
                                     auto integrator,
                                     const std::vector<double>& dts) {
        std::cout << "------------------------------------------------------------\n";
        std::cout << name << " (theoretical order " << theoretical_order << ")\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << std::left
                  << std::setw(12) << "dt (s)"
                  << std::setw(15) << "Error (m)"
                  << std::setw(15) << "Error ratio"
                  << "Emp. order\n";

        double prev_error = -1.0;
        for (double dt : dts) {
            const State4D result = run_to_t_end(integrator, t_end, dt, s0, deriv);
            const double err = position_error(result, truth);

            std::cout << std::left
                      << std::defaultfloat << std::setprecision(6)
                      << std::setw(12) << dt
                      << std::scientific << std::setprecision(3)
                      << std::setw(15) << err;

            if (prev_error > 0.0 && err > 0.0) {
                const double ratio = prev_error / err;
                const double order = std::log2(ratio);
                std::cout << std::defaultfloat << std::setprecision(3)
                          << std::setw(15) << ratio
                          << std::fixed << std::setprecision(2) << order;
            } else {
                std::cout << std::setw(15) << "-" << "-";
            }
            std::cout << "\n";
            prev_error = err;
        }
        std::cout << "\n";
    };

    run_convergence_table("Euler",        1, euler_step<State4D, decltype(deriv)>, dts_low_order);
    run_convergence_table("Heun",         2, heun_step<State4D, decltype(deriv)>,  dts_low_order);
    run_convergence_table("RK4",          4, rk4_step<State4D, decltype(deriv)>,   dts_rk4);
    run_convergence_table("RK45 (fixed)", 5, rk45_step<State4D, decltype(deriv)>,  dts_rk45);
    run_convergence_table("RK8",          8, rk8_step<State4D, decltype(deriv)>,   dts_rk8);

    return 0;
}
