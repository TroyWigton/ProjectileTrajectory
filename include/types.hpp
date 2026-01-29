/**
 * @file types.hpp
 * @brief Defines common data types and indices used throughout the simulation.
 */

#ifndef TYPES_HPP
#define TYPES_HPP

#include <array>
#include <functional>

/**
 * @namespace StateIndex4D
 * @brief Indices for accessing components of a 4D state vector (classic projectile).
 */
namespace StateIndex4D {
    /** Enumeration for 4D state vector indices */
    enum { 
        X_POS = 0, ///< X Position index
        Y_POS = 1, ///< Y Position index
        X_VEL = 2, ///< X Velocity index
        Y_VEL = 3  ///< Y Velocity index
    };
}

/**
 * @namespace StateIndex6D
 * @brief Indices for accessing components of a 6D state vector (3D motion).
 */
namespace StateIndex6D {
    /** Enumeration for 6D state vector indices */
    enum { 
        X_POS = 0, ///< X Position index
        Y_POS = 1, ///< Y Position index
        Z_POS = 2, ///< Z Position index
        X_VEL = 3, ///< X Velocity index
        Y_VEL = 4, ///< Y Velocity index
        Z_VEL = 5  ///< Z Velocity index
    };
}

/**
 * @brief Represents a 4D state vector [x_pos, y_pos, x_vel, y_vel].
 */
using State4D = std::array<double, 4>;

/**
 * @brief Represents a 6D state vector [x, y, z, vx, vy, vz].
 */
using State6D = std::array<double, 6>;

/**
 * @brief Function pointer type for calculating system derivatives.
 * 
 * Defines a specific physics model (forces/derivatives) including physical constants. 
 * This is bound into a SystemDerivative function within the Simulation class.
 */
using DerivativeFuncPtr = std::function<void(const State4D& current_state, double current_time, State4D& derivative_output, double gravity, double k_over_m)>;


#endif // TYPES_HPP