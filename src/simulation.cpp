#include "simulation.hpp"
#include "constants.hpp"
#include <cmath>

Simulation::Simulation(double v0,
                                             DerivativeFuncPtr derivative,
                                             IntegratorWithParamsFuncPtr integrator,
                                             double g,
                                             double k_over_m,
                                             double h0)
        : Simulation(g, k_over_m, v0, h0, derivative, integrator) {}

Simulation::Simulation(double g,
                                             double k_over_m,
                                             double v0,
                                             double h0,
                                             DerivativeFuncPtr derivative,
                                             IntegratorWithParamsFuncPtr integrator)
        : m_g(g),
            m_k_over_m(k_over_m),
            m_v0(v0),
            m_h0(h0),
            m_derivative(derivative),
            m_integrator(integrator) {}

ScenarioResult Simulation::run(double angle_deg, double dt) const {
    const double angle_rad = angle_deg * M_PI / 180.0;
    State state = {0.0, m_h0, m_v0 * std::cos(angle_rad), m_v0 * std::sin(angle_rad)};
    double t = 0.0;
    double last_x = 0.0;
    double last_y = m_h0;

    while (state[Y_POS] >= 0.0) {
        last_x = state[X_POS];
        last_y = state[Y_POS];
        state = m_integrator(state, t, dt, m_derivative, m_g, m_k_over_m);
        t += dt;
    }

    const double denominator = state[Y_POS] - last_y;
    const double impact_x = std::fabs(denominator) > 1e-12
                                ? last_x + (-last_y) * (state[X_POS] - last_x) / denominator
                                : state[X_POS];

    return {angle_deg, m_v0, m_h0, dt, impact_x};
}
