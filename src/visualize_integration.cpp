/**
 * @file visualize_integration.cpp
 * @brief Generates CSV datasets that visualize numerical integration behavior.
 *
 * Produces the raw data behind the integration_visualization.ipynb notebook,
 * answering three questions for the canonical golf-ball v^2-drag trajectory:
 *
 *   1. Step size:  how does a single integrator's accuracy improve as dt shrinks?
 *   2. Method:     how do different integrators compare at the *same* dt?
 *   3. Stages:     how do the "minor step" derivative probes of a multi-stage
 *                  Runge-Kutta method combine into a "major step" estimate that
 *                  is more accurate than any single slope -- shrinking the error
 *                  for a fixed step size as the method order rises?
 *
 * The trajectory, error, and convergence studies use the canonical
 * drag_deriv_v_squared model (derivative_functions.hpp); the stage-decomposition
 * study uses a 1-D sample ODE (y' = y) so each stage slope is a single number.
 * Both drive the project's existing integrators (integrators.hpp). The stage
 * logging re-expresses the Euler/Heun/RK4 stage structure explicitly so every
 * intermediate slope can be recorded; the reconstructed end-of-step state is
 * asserted equal to the canonical *_step() output, guaranteeing the logged
 * stages faithfully represent the real integrators.
 *
 * Output (written to argv[1], default "viz_data"):
 *   meta.csv          scenario constants (key,value)
 *   trajectories.csv  full arcs: step-size study, integrator study, reference
 *   error_vs_time.csv position error vs time, per integrator, at a fixed dt
 *   convergence.csv   global error at fixed T vs dt, per integrator
 *   scalar_stages.csv per-stage eval points + slopes for one step of y' = y
 *   scalar_summary.csv per-method endpoint, blended slope, secant slope, error
 *   scalar_truth.csv  dense true solution y = e^t over the step window
 */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <cstdlib>
#include "../include/constants.hpp"
#include "../include/types.hpp"
#include "../include/integrators.hpp"
#include "../include/derivative_functions.hpp"

using namespace StateIndex4D;

struct TrajPoint { double t, x, y, vx, vy; };

// --- generic stepping helpers ----------------------------------------------

// Integrate from s0 with fixed dt until the projectile crosses y < 0 (or a
// safety cap). Records the state after every step, including the first point
// below ground so the plotted arc closes cleanly at the axis.
template<typename Integrator, typename Deriv>
std::vector<TrajPoint> integrate_to_ground(Integrator integ, double dt,
                                           const State4D& s0, Deriv deriv,
                                           double t_cap = 40.0) {
    std::vector<TrajPoint> pts;
    State4D s = s0;
    double t = 0.0;
    pts.push_back({t, s[X_POS], s[Y_POS], s[X_VEL], s[Y_VEL]});
    while (s[Y_POS] >= 0.0 && t < t_cap) {
        s = integ(s, t, dt, deriv);
        t += dt;
        pts.push_back({t, s[X_POS], s[Y_POS], s[X_VEL], s[Y_VEL]});
    }
    return pts;
}

// Dense uniform-dt integration to a fixed T (recomputing dt to land on T),
// recording every step. Used to build the high-precision reference arc.
template<typename Integrator, typename Deriv>
std::vector<TrajPoint> dense_to_T(Integrator integ, double T, double dt,
                                  const State4D& s0, Deriv deriv,
                                  double& actual_dt_out) {
    const long long steps = std::lround(T / dt);
    const double adt = T / static_cast<double>(steps);
    actual_dt_out = adt;
    std::vector<TrajPoint> pts;
    pts.reserve(static_cast<size_t>(steps) + 1);
    State4D s = s0;
    double t = 0.0;
    pts.push_back({0.0, s[X_POS], s[Y_POS], s[X_VEL], s[Y_VEL]});
    for (long long i = 0; i < steps; ++i) {
        s = integ(s, t, adt, deriv);
        t += adt;
        pts.push_back({t, s[X_POS], s[Y_POS], s[X_VEL], s[Y_VEL]});
    }
    return pts;
}

// Linear interpolation of a uniform-dt reference trajectory at arbitrary t.
TrajPoint ref_at(const std::vector<TrajPoint>& ref, double adt, double t) {
    if (t <= 0.0) return ref.front();
    const double fidx = t / adt;
    const size_t i = static_cast<size_t>(std::floor(fidx));
    if (i + 1 >= ref.size()) return ref.back();
    const double f = fidx - static_cast<double>(i);
    const TrajPoint& a = ref[i];
    const TrajPoint& b = ref[i + 1];
    return { t,
             a.x  + f * (b.x  - a.x),
             a.y  + f * (b.y  - a.y),
             a.vx + f * (b.vx - a.vx),
             a.vy + f * (b.vy - a.vy) };
}

// Integrate to a fixed T and return only the final state (convergence study).
template<typename Integrator, typename Deriv>
State4D run_to_T(Integrator integ, double T, double dt, const State4D& s0, Deriv deriv) {
    const long long steps = std::lround(T / dt);
    const double adt = T / static_cast<double>(steps);
    State4D s = s0;
    double t = 0.0;
    for (long long i = 0; i < steps; ++i) { s = integ(s, t, adt, deriv); t += adt; }
    return s;
}

double pos_err(const State4D& a, const State4D& b) {
    const double dx = a[X_POS] - b[X_POS], dy = a[Y_POS] - b[Y_POS];
    return std::sqrt(dx * dx + dy * dy);
}

// --- scalar stage decomposition (minor-step slopes -> major-step slope) ------
// To show the idea cleanly we use a 1-D sample ODE, y' = f(y), where each stage
// slope is a single number. Each ScalarStage records the point (eval_t, eval_y)
// where the derivative was sampled, the slope found there, and this stage's
// weight in the final blend. The reconstructed endpoint is asserted equal to the
// canonical *_step() output (run on a FixedState<1>), so these are exactly the
// stages the real integrators take.

struct ScalarStage {
    double eval_t;   // time at which the derivative was sampled
    double eval_y;   // predicted y at which the derivative was sampled
    double slope;    // the derivative (dy/dt) found there
    double weight;   // weight of this stage in the final combination
};
struct ScalarRun { std::vector<ScalarStage> stages; double y1; };

template<typename F>
ScalarRun euler_scalar(double t0, double y0, double h, F f) {
    const double k1 = f(y0);
    return { { {t0, y0, k1, 1.0} }, y0 + h * k1 };
}

template<typename F>
ScalarRun heun_scalar(double t0, double y0, double h, F f) {
    const double k1 = f(y0);
    const double yp = y0 + h * k1;            // Euler predictor to the step end
    const double k2 = f(yp);
    return { { {t0, y0, k1, 0.5}, {t0 + h, yp, k2, 0.5} },
             y0 + 0.5 * h * (k1 + k2) };
}

template<typename F>
ScalarRun rk4_scalar(double t0, double y0, double h, F f) {
    const double k1 = f(y0);                  // slope at the start (major point)
    const double ya = y0 + 0.5 * h * k1;
    const double k2 = f(ya);                  // slope at the midpoint (minor)
    const double yb = y0 + 0.5 * h * k2;
    const double k3 = f(yb);                  // refined midpoint slope (minor)
    const double yc = y0 + h * k3;
    const double k4 = f(yc);                  // slope at the end (major point)
    const double y1 = y0 + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    return { { {t0,           y0, k1, 1.0 / 6.0},
               {t0 + 0.5 * h, ya, k2, 2.0 / 6.0},
               {t0 + 0.5 * h, yb, k3, 2.0 / 6.0},
               {t0 + h,       yc, k4, 1.0 / 6.0} }, y1 };
}

// --- CSV utilities ----------------------------------------------------------
std::ofstream open_csv(const std::string& dir, const std::string& name) {
    std::ofstream f(dir + "/" + name);
    if (!f.is_open()) {
        std::cerr << "ERROR: cannot open " << dir << "/" << name
                  << " (does the output directory exist/writable?)\n";
        std::exit(1);
    }
    f << std::setprecision(12);
    return f;
}

int main(int argc, char** argv) {
    const std::string outdir = (argc > 1) ? argv[1] : "viz_data";
    // The output directory must already exist; open_csv() reports clearly if not.

    // --- scenario: a lofted v^2-drag shot --------------------------------
    // Chosen to be smooth and long-lived: moderate launch speed and drag
    // (a_drag ~ a few m/s^2, gentler than gravity) so coarse steps stay
    // well-behaved and error accumulates gradually over a ~9 s flight, rather
    // than dumping almost all of it into a stiff first step.
    const double g          = GRAVITY_EARTH;
    const double k_over_m    = 0.0025;   // moderate drag (vs 0.0057 golf at v0=100)
    const double v0          = 60.0;     // slower launch -> smoother, longer flight
    const double launch_deg  = 60.0;     // high arc -> long time aloft, speeds up on descent
    const double launch_rad  = launch_deg * M_PI / 180.0;
    const State4D s0 = {0.0, 0.0, v0 * std::cos(launch_rad), v0 * std::sin(launch_rad)};

    ProjectileContext ctx{g, k_over_m};
    auto deriv = [&ctx](const State4D& s, double t, State4D& d) {
        drag_deriv_v_squared(s, t, d, ctx);
    };

    // Existing project integrators, instantiated for this state+deriv type.
    auto E   = euler_step<State4D, decltype(deriv)>;
    auto H   = heun_step<State4D, decltype(deriv)>;
    auto R4  = rk4_step<State4D, decltype(deriv)>;
    auto R45 = rk45_step<State4D, decltype(deriv)>;
    auto R8  = rk8_step<State4D, decltype(deriv)>;

    // --- dense reference arc (truth) -------------------------------------
    const double T_ref  = 16.0;     // safely beyond the drag-shortened flight
    const double ref_dt = 1.0e-4;
    double ref_adt = 0.0;
    const std::vector<TrajPoint> ref = dense_to_T(R8, T_ref, ref_dt, s0, deriv, ref_adt);

    // ===== meta.csv =====
    {
        std::ofstream f = open_csv(outdir, "meta.csv");
        f << "key,value\n"
          << "model,v_squared_drag\n"
          << "g," << g << "\n"
          << "k_over_m," << k_over_m << "\n"
          << "v0," << v0 << "\n"
          << "launch_deg," << launch_deg << "\n"
          << "ref_integrator,RK8\n"
          << "ref_dt," << ref_adt << "\n";
    }

    // ===== trajectories.csv =====
    const double dt_traj = 1.5;  // coarse dt for the integrator-comparison trajectory + zoom
    {
        std::ofstream f = open_csv(outdir, "trajectories.csv");
        f << "study,integrator,dt,step,t,x,y,vx,vy\n";
        auto dump = [&](const std::string& study, const std::string& integ,
                        double dt, const std::vector<TrajPoint>& pts) {
            for (size_t i = 0; i < pts.size(); ++i)
                f << study << "," << integ << "," << dt << "," << i << ","
                  << pts[i].t << "," << pts[i].x << "," << pts[i].y << ","
                  << pts[i].vx << "," << pts[i].vy << "\n";
        };

        // (1) Step-size study: one method (Euler), several dt -> see it converge.
        const double stepsize_dts[] = {1.0, 0.5, 0.25, 0.1};
        for (double dt : stepsize_dts)
            dump("stepsize", "Euler", dt, integrate_to_ground(E, dt, s0, deriv));

        // (2) Integrator study: all methods at one coarse dt -> compare orders.
        dump("integrator", "Euler", dt_traj, integrate_to_ground(E,   dt_traj, s0, deriv));
        dump("integrator", "Heun",  dt_traj, integrate_to_ground(H,   dt_traj, s0, deriv));
        dump("integrator", "RK4",   dt_traj, integrate_to_ground(R4,  dt_traj, s0, deriv));
        dump("integrator", "RK45",  dt_traj, integrate_to_ground(R45, dt_traj, s0, deriv));
        dump("integrator", "RK8",   dt_traj, integrate_to_ground(R8,  dt_traj, s0, deriv));

        // (3) Reference arc: downsampled dense RK8, stopping just past ground.
        const int stride = std::max(1, static_cast<int>(std::lround(0.02 / ref_adt)));
        for (size_t i = 0; i < ref.size(); i += static_cast<size_t>(stride)) {
            f << "reference,RK8_ref," << ref_adt << "," << i << ","
              << ref[i].t << "," << ref[i].x << "," << ref[i].y << ","
              << ref[i].vx << "," << ref[i].vy << "\n";
            if (ref[i].y < 0.0) break;
        }
    }

    // ===== error_vs_time.csv (all integrators at a shared moderate dt) =====
    {
        const double dt_err = 0.1;  // fine enough that error builds gradually over the flight (not first-step-dominated)
        std::ofstream f = open_csv(outdir, "error_vs_time.csv");
        f << "integrator,dt,t,x,y,x_ref,y_ref,abs_error\n";
        auto dump_err = [&](const std::string& name, const std::vector<TrajPoint>& pts) {
            for (const auto& p : pts) {
                if (p.t > T_ref) break;
                const TrajPoint r = ref_at(ref, ref_adt, p.t);
                const double e = std::sqrt((p.x - r.x) * (p.x - r.x) +
                                           (p.y - r.y) * (p.y - r.y));
                f << name << "," << dt_err << "," << p.t << "," << p.x << "," << p.y
                  << "," << r.x << "," << r.y << "," << e << "\n";
            }
        };
        dump_err("Euler", integrate_to_ground(E,   dt_err, s0, deriv));
        dump_err("Heun",  integrate_to_ground(H,   dt_err, s0, deriv));
        dump_err("RK4",   integrate_to_ground(R4,  dt_err, s0, deriv));
        dump_err("RK45",  integrate_to_ground(R45, dt_err, s0, deriv));
        dump_err("RK8",   integrate_to_ground(R8,  dt_err, s0, deriv));
    }

    // ===== convergence.csv (global error at fixed T vs dt) =====
    // The convergence study uses a more demanding shot than the display scenario
    // above (faster launch, stronger drag). The smooth lofted shot is so easy
    // that the high-order methods hit the floating-point round-off floor almost
    // immediately, hiding their order; this faster problem keeps RK4/RK45/RK8 in
    // their asymptotic regime across the dt sweep so the slopes stay observable.
    {
        ProjectileContext ctx_c{GRAVITY_EARTH, 0.0057};
        auto deriv_c = [&ctx_c](const State4D& s, double t, State4D& d) {
            drag_deriv_v_squared(s, t, d, ctx_c);
        };
        const double v0c = 100.0, angc = 45.0 * M_PI / 180.0;
        const State4D s0c = {0.0, 0.0, v0c * std::cos(angc), v0c * std::sin(angc)};
        auto Ec   = euler_step<State4D, decltype(deriv_c)>;
        auto Hc   = heun_step<State4D, decltype(deriv_c)>;
        auto R4c  = rk4_step<State4D, decltype(deriv_c)>;
        auto R45c = rk45_step<State4D, decltype(deriv_c)>;
        auto R8c  = rk8_step<State4D, decltype(deriv_c)>;

        const double T = 5.0;
        const double conv_ref_dt = 1.0e-5;  // RK8 reference; ~round-off-floor accurate
        const State4D truth = run_to_T(R8c, T, conv_ref_dt, s0c, deriv_c);

        std::ofstream f = open_csv(outdir, "convergence.csv");
        f << "integrator,theoretical_order,dt,error\n";
        auto sweep = [&](const std::string& name, int order, auto integ,
                         const std::vector<double>& dts) {
            for (double dt : dts) {
                const State4D r = run_to_T(integ, T, dt, s0c, deriv_c);
                f << name << "," << order << "," << dt << "," << pos_err(r, truth) << "\n";
            }
        };
        // One common, wide dt range for every method, so the lines are directly
        // comparable and any large-dt crossing of the steep (high-order) lines
        // over the shallow (low-order) ones would be visible.
        const std::vector<double> common = {2.0, 1.0, 0.5, 0.25, 0.125,
                                             0.0625, 0.03125, 0.015625};
        sweep("Euler", 1, Ec,   common);
        sweep("Heun",  2, Hc,   common);
        sweep("RK4",   4, R4c,  common);
        sweep("RK45",  5, R45c, common);
        sweep("RK8",   8, R8c,  common);
    }

    // ===== scalar_stages.csv + scalar_summary.csv + scalar_truth.csv =====
    // One step of the sample ODE y' = y from y0 = 1 (true solution y = e^t),
    // dissected so the stage slopes and their weighted blend are visible. The
    // step is run through the project's real integrators on a FixedState<1>;
    // the explicitly-logged stages are asserted to reproduce that output.
    {
        using Scalar = FixedState<1>;
        auto sderiv  = [](const Scalar& s, double, Scalar& d) { d[0] = s[0]; }; // y' = y
        auto fscalar = [](double y) { return y; };
        auto SE  = euler_step<Scalar, decltype(sderiv)>;
        auto SH  = heun_step<Scalar, decltype(sderiv)>;
        auto SR4 = rk4_step<Scalar, decltype(sderiv)>;

        const double t0 = 0.0, y0 = 1.0, h = 1.0;
        const double t1 = t0 + h;
        const double y_true = y0 * std::exp(t1);       // analytic truth for y' = y

        std::ofstream fs = open_csv(outdir, "scalar_stages.csv");
        fs << "method,h,stage,eval_t,eval_y,slope,weight\n";
        std::ofstream fsum = open_csv(outdir, "scalar_summary.csv");
        fsum << "method,h,t0,y0,t1,y1,k_eff,y_true,secant_slope,error\n";

        auto emit = [&](const std::string& name, const ScalarRun& r, double canon_y1) {
            if (std::abs(r.y1 - canon_y1) > 1e-12) {
                std::cerr << "FIDELITY FAIL: " << name << " reconstructed endpoint "
                          << r.y1 << " vs canonical " << canon_y1 << "\n";
                std::exit(2);
            }
            for (size_t i = 0; i < r.stages.size(); ++i) {
                const ScalarStage& st = r.stages[i];
                fs << name << "," << h << "," << i << "," << st.eval_t << ","
                   << st.eval_y << "," << st.slope << "," << st.weight << "\n";
            }
            const double k_eff  = (r.y1 - y0) / h;          // the major-step slope
            const double secant = (y_true - y0) / h;        // exact average slope
            fsum << name << "," << h << "," << t0 << "," << y0 << "," << t1 << ","
                 << r.y1 << "," << k_eff << "," << y_true << "," << secant << ","
                 << std::abs(r.y1 - y_true) << "\n";
        };

        const Scalar ss0 = {y0};
        emit("Euler", euler_scalar(t0, y0, h, fscalar), SE(ss0, t0, h, sderiv)[0]);
        emit("Heun",  heun_scalar(t0, y0, h, fscalar),  SH(ss0, t0, h, sderiv)[0]);
        emit("RK4",   rk4_scalar(t0, y0, h, fscalar),   SR4(ss0, t0, h, sderiv)[0]);

        // Dense true solution over a window around the step, for the plot.
        std::ofstream ft = open_csv(outdir, "scalar_truth.csv");
        ft << "t,y\n";
        const double ta = t0 - 0.12 * h, tb = t1 + 0.28 * h;
        const int N = 300;
        for (int j = 0; j <= N; ++j) {
            const double tt = ta + (tb - ta) * (static_cast<double>(j) / N);
            ft << tt << "," << y0 * std::exp(tt) << "\n";
        }
    }

    std::cout << "Wrote CSVs to '" << outdir << "/': meta, trajectories, error_vs_time, "
                 "convergence, scalar_stages, scalar_summary, scalar_truth\n";
    return 0;
}
