#ifndef SIMULATION_HPP
#define SIMULATION_HPP

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

class Simulation {
public:
    Simulation(double k_over_m,
               SystemIntegrator integrator,
               double v0 = 100.0,
               DerivativeFuncPtr derivative = drag_deriv_v_squared,
               double g = 9.81,
               double h0 = 0.0);

    ScenarioResult run(double angle_deg, double dt) const;

    double gravity() const { return g_; }
    double drag_ratio() const { return k_over_m_; }
    double initial_velocity() const { return v0_; }
    double initial_height() const { return h0_; }

private:
    double g_;
    double k_over_m_;
    double v0_;
    double h0_;
    
    // Abstracted stepper function that knows how to advance the system
    // (Integrator + Derivative + Parameters all bound)
    std::function<State(const State&, double, double)> stepper_;
};

#endif // SIMULATION_HPP
