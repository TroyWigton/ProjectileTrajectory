#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <cmath>
#include "types.hpp"
#include "integrators.hpp"
#include "derivative_functions.hpp"

struct ScenarioResult {
    double angle;
    double velocity;
    double height;
    double deltaT;
    double distance;
    long long step_count;
};

// Binds one integration method (Integrator) to one physics model
// (ProjectileDerivative) and a fixed set of physical parameters, then owns a
// single projectile's state so it can be stepped forward in time.
//
// Integrator is a compile-time template parameter — typically a stateless
// functor with a templated call operator (EulerStepper, RK4Stepper, ... from
// integrators.hpp) — rather than a std::function. The derivative binding
// inside step() is a concrete lambda type resolved once per
// Simulation<Integrator> instantiation, so nothing is boxed on the per-step
// hot path and the integrator+derivative call chain inlines.
//
// Integrator can still be a genuine std::function (SystemIntegrator<State4D>)
// when a caller actually needs to pick the method at runtime — e.g. looping
// over several integrators to compare them (see test.cpp) — that case just
// isn't paid for by callers who don't need it.
template <typename Integrator>
class Simulation {
public:
    Simulation(double k_over_m,
               Integrator integrator,
               double v0 = 100.0,
               ProjectileDerivative derivative = drag_deriv_v_squared,
               double g = 9.81,
               double h0 = 0.0)
        : g_(g), k_over_m_(k_over_m), v0_(v0), h0_(h0),
          ctx_{g, k_over_m}, integrator_(integrator), derivative_(derivative),
          state_{0.0, h0, 0.0, 0.0}, t_(0.0) {}

    // Advances the owned state by one step of size dt.
    void step(double dt) {
        auto bound_deriv = [this](const State4D& s, double time, State4D& out) {
            derivative_(s, time, out, ctx_);
        };
        state_ = integrator_(state_, t_, dt, bound_deriv);
        t_ += dt;
    }

    // Launches at angle_deg and steps to ground impact (y < 0), linearly
    // interpolating between the last two states for the exact impact point.
    ScenarioResult run(double angle_deg, double dt) {
        using namespace StateIndex4D;
        const double angle_rad = angle_deg * M_PI / 180.0;
        state_ = {0.0, h0_, v0_ * std::cos(angle_rad), v0_ * std::sin(angle_rad)};
        t_ = 0.0;
        double last_x = 0.0;
        double last_y = h0_;
        long long steps = 0;

        while (state_[Y_POS] >= 0.0) {
            last_x = state_[X_POS];
            last_y = state_[Y_POS];
            step(dt);
            steps++;
        }

        const double denominator = state_[Y_POS] - last_y;
        const double impact_x = std::fabs(denominator) > 1e-12
                                    ? last_x + (-last_y) * (state_[X_POS] - last_x) / denominator
                                    : state_[X_POS];

        return {angle_deg, v0_, h0_, dt, impact_x, steps};
    }

    const State4D& state() const { return state_; }
    double time() const { return t_; }

    double gravity() const { return g_; }
    double drag_ratio() const { return k_over_m_; }
    double initial_velocity() const { return v0_; }
    double initial_height() const { return h0_; }

private:
    const double g_;
    const double k_over_m_;
    const double v0_;
    const double h0_;
    const ProjectileContext ctx_;
    Integrator integrator_;
    ProjectileDerivative derivative_;
    State4D state_;
    double t_;
};

// Deduces Integrator so callers don't have to spell out the template
// argument (class templates can't deduce from constructor args until C++17,
// and this project targets C++14): auto sim = make_simulation(k_over_m, RK4Stepper{}, v0);
template <typename Integrator>
Simulation<Integrator> make_simulation(double k_over_m,
                                        Integrator integrator,
                                        double v0 = 100.0,
                                        ProjectileDerivative derivative = drag_deriv_v_squared,
                                        double g = 9.81,
                                        double h0 = 0.0) {
    return Simulation<Integrator>(k_over_m, integrator, v0, derivative, g, h0);
}

#endif // SIMULATION_HPP
