#ifndef DERIVATIVE_FUNCTIONS_HPP
#define DERIVATIVE_FUNCTIONS_HPP

#include "types.hpp"

void drag_deriv_v_squared(const State& s, double t, State& deriv, double g, double k_over_m);

void drag_deriv_linear(const State& s, double t, State& deriv, double g, double k_over_m);

void no_drag_deriv(const State& s, double t, State& deriv, double g, double k_over_m = 0.0);

#endif // DERIVATIVE_FUNCTIONS_HPP
