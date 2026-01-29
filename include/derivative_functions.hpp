/**
 * @file derivative_functions.hpp
 * @brief Declarations of derivative functions (physics models) for the simulation.
 */

#ifndef DERIVATIVE_FUNCTIONS_HPP
#define DERIVATIVE_FUNCTIONS_HPP

#include "types.hpp"

/**
 * @brief Computes derivatives for quadratic drag (v^2).
 * 
 * Modeled as F_drag = -k * v^2.
 * 
 * @param s Current state vector.
 * @param t Current time.
 * @param deriv Output state vector for derivatives.
 * @param g Gravitational acceleration.
 * @param k_over_m Drag coefficient divided by mass.
 */
void drag_deriv_v_squared(const State4D& s, double t, State4D& deriv, double g, double k_over_m);

/**
 * @brief Computes derivatives for linear drag (v).
 * 
 * Modeled as F_drag = -k * v.
 * 
 * @param s Current state vector.
 * @param t Current time.
 * @param deriv Output state vector for derivatives.
 * @param g Gravitational acceleration.
 * @param k_over_m Drag coefficient divided by mass.
 */
void drag_deriv_linear(const State4D& s, double t, State4D& deriv, double g, double k_over_m);

/**
 * @brief Computes derivatives for no drag (vacuum).
 * 
 * Only gravity acts on the projectile.
 * 
 * @param s Current state vector.
 * @param t Current time.
 * @param deriv Output state vector for derivatives.
 * @param g Gravitational acceleration.
 * @param k_over_m Unused (defaults to 0).
 */
void no_drag_deriv(const State4D& s, double t, State4D& deriv, double g, double k_over_m = 0.0);

/**
 * @brief Computes derivatives for golf ball physics.
 * 
 * Uses variable drag coefficient based on Reynolds number / velocity.
 * 
 * @param s Current state vector.
 * @param t Current time.
 * @param deriv Output state vector for derivatives.
 * @param g Gravitational acceleration.
 * @param v_critical Critical velocity parameter (unused/legacy in some impls).
 */
void golf_ball_drag_deriv(const State4D& s, double t, State4D& deriv, double g, double v_critical);

/**
 * @brief Computes derivatives for a variable drag model.
 * 
 * @param s Current state vector.
 * @param t Current time.
 * @param deriv Output state vector for derivatives.
 * @param g Gravitational acceleration.
 * @param base_k_over_m Base drag coefficient factor.
 */
void variable_drag_deriv(const State4D& s, double t, State4D& deriv, double g, double base_k_over_m);

#endif // DERIVATIVE_FUNCTIONS_HPP
