#ifndef TYPES_HPP
#define TYPES_HPP

#include <array>
#include <functional>

// 4D State Indices (Namespaced)
namespace StateIndex4D {
    enum { X_POS = 0, Y_POS = 1, X_VEL = 2, Y_VEL = 3 };
}

// 6D State Indices (Namespaced)
namespace StateIndex6D {
    enum { X_POS = 0, Y_POS = 1, Z_POS = 2, X_VEL = 3, Y_VEL = 4, Z_VEL = 5 };
}

// 4D state: [x_pos, y_pos, x_vel, y_vel]
using State4D = std::array<double, 4>;

// 6D state: [x, y, z, vx, vy, vz]
using State6D = std::array<double, 6>;

// Specific physics model (forces/derivatives) including physical constants. 
// Bound into a SystemDerivative within the Simulation class.
using DerivativeFuncPtr = std::function<void(const State4D& current_state, double current_time, State4D& derivative_output, double gravity, double k_over_m)>;

#endif // TYPES_HPP