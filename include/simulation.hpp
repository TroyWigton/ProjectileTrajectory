#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "types.hpp"
#include "integrators.hpp"

struct ScenarioResult {
    double angle;
    double velocity;
    double height;
    double deltaT;
    double distance;
};

class Simulation {
public:
    Simulation(double v0,
               DerivativeFuncPtr derivative,
               IntegratorWithParamsFuncPtr integrator = rk4_step_with_params,
               double g = 9.81,
               double k_over_m = 0.0,
               double h0 = 0.0);

    Simulation(double g,
               double k_over_m,
               double v0,
               double h0,
               DerivativeFuncPtr derivative,
               IntegratorWithParamsFuncPtr integrator = rk4_step_with_params);

    ScenarioResult run(double angle_deg, double dt) const;

    double gravity() const { return m_g; }
    double drag_ratio() const { return m_k_over_m; }
    double initial_velocity() const { return m_v0; }
    double initial_height() const { return m_h0; }

private:
    double m_g;
    double m_k_over_m;
    double m_v0;
    double m_h0;
    DerivativeFuncPtr m_derivative;
    IntegratorWithParamsFuncPtr m_integrator;
};

#endif // SIMULATION_HPP
