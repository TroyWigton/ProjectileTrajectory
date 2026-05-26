#ifndef DERIVATIVE_FUNCTIONS_HPP
#define DERIVATIVE_FUNCTIONS_HPP

#include "types.hpp"

// Problem-specific context bundle for projectile-motion derivatives.
// Bundled into the generic ParameterizedDerivative<State, Context> signature
// so the framework never sees projectile-specific parameters.
struct ProjectileContext {
    double g;        // gravitational acceleration (m/s^2)
    double k_over_m; // drag coefficient ratio (k/m)
};

// Convenience alias for the projectile-motion derivative signature.
using ProjectileDerivative = ParameterizedDerivative<State4D, ProjectileContext>;

void drag_deriv_v_squared(const State4D& s, double t, State4D& deriv, const ProjectileContext& ctx);

void drag_deriv_linear(const State4D& s, double t, State4D& deriv, const ProjectileContext& ctx);

void no_drag_deriv(const State4D& s, double t, State4D& deriv, const ProjectileContext& ctx);

void golf_ball_drag_deriv(const State4D& s, double t, State4D& deriv, const ProjectileContext& ctx);

void variable_drag_deriv(const State4D& s, double t, State4D& deriv, const ProjectileContext& ctx);

#endif // DERIVATIVE_FUNCTIONS_HPP
