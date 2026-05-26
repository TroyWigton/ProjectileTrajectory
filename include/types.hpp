#ifndef TYPES_HPP
#define TYPES_HPP

#include <array>
#include <vector>
#include <functional>

// State containers ----------------------------------------------------------
// The framework treats the simulation state as an opaque flat container of
// doubles; interpretation of the slots is owned by each problem domain.
//
// FixedState<N>: compile-time size, stack-allocated, trivially copyable.
//                Use when N is small and known at build time (projectile
//                motion, demos, anything where N is a handful).
// DynamicState:  runtime-sized, heap-backed. Use when N is determined at
//                runtime, varies between runs of the same binary, or is too
//                large to live on the stack.
//
// Integrator templates iterate via .size() and operator[], so they accept
// either container uniformly.

template <std::size_t N>
using FixedState = std::array<double, N>;

using DynamicState = std::vector<double>;

// Projectile-motion state shapes used in this repo. New problem domains
// should declare their own typedef in terms of FixedState<N> or DynamicState.
using State4D = FixedState<4>;  // [x, y, vx, vy]
using State6D = FixedState<6>;  // [x, y, z, vx, vy, vz]


// State interpretation (problem-specific) -----------------------------------
// Anonymous-enum-in-namespace gives named, scoped integer indices into the
// flat state array. Each problem domain provides its own mapping; the
// framework itself never reads these.

namespace StateIndex4D {
    enum { X_POS = 0, Y_POS = 1, X_VEL = 2, Y_VEL = 3 };
}

namespace StateIndex6D {
    enum { X_POS = 0, Y_POS = 1, Z_POS = 2, X_VEL = 3, Y_VEL = 4, Z_VEL = 5 };
}


// Derivative-function signatures --------------------------------------------
// Generic, interpretation-neutral form: takes (state, t, deriv_out, ctx)
// where Context is whatever problem-specific parameter bundle the derivative
// needs (gravity + drag for projectile motion; constants + masses for an
// N-body system; anything else for some other problem). The framework never
// inspects Context — it just passes it through to the bound derivative.

template <typename State, typename Context>
using ParameterizedDerivative =
    std::function<void(const State& state, double t, State& deriv_out, const Context& ctx)>;

#endif // TYPES_HPP
