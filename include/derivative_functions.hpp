#ifndef DERIVATIVE_FUNCTIONS_HPP
#define DERIVATIVE_FUNCTIONS_HPP

#include "types.hpp"

void drag_deriv_v_squared(const State4D& s, double t, State4D& deriv, double g, double k_over_m);

void drag_deriv_linear(const State4D& s, double t, State4D& deriv, double g, double k_over_m);

void no_drag_deriv(const State4D& s, double t, State4D& deriv, double g, double k_over_m = 0.0);

void golf_ball_drag_deriv(const State4D& s, double t, State4D& deriv, double g, double);

void variable_drag_deriv(const State4D& s, double t, State4D& deriv, double g, double base_k_over_m);

#endif // DERIVATIVE_FUNCTIONS_HPP
