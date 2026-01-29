/**
 * @file constants.hpp
 * @brief Defines physical constants and simulation parameters.
 */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/// Acceleration due to gravity on Earth (m/s^2)
constexpr double GRAVITY_EARTH = 9.81;

/**
 * @namespace DragRatios
 * @brief Pre-defined drag-to-mass ratios (k/m) for various projectiles.
 * 
 * These values approximate the ballistic coefficient factors used in the simulation.
 */
namespace DragRatios {
    constexpr double VACUUM = 0.0;          ///< No drag
    constexpr double BULLET_223 = 0.0012;   ///< Approx for .223 Rem (55gr), derived from Cd ~ 0.23
    constexpr double BULLET_9MM = 0.0015;   ///< Approx for 9mm projectile
    constexpr double GOLF_BALL = 0.0057;    ///< Standard golf ball approximation
    constexpr double PING_PONG_BALL = 0.134;///< Ping Pong ball (high drag)
}

/**
 * @namespace GolfBallPhysics
 * @brief Specific physical parameters for detailed golf ball simulation.
 * 
 * Includes mass, diameter, and aerodynamic properties for variable drag models.
 */
namespace GolfBallPhysics {
    // Physical properties
    constexpr double MASS_KG = 0.04593;     ///< Mass of a standard golf ball (kg)
    constexpr double DIAMETER_M = 0.04267;  ///< Diameter of a standard golf ball (m)
    constexpr double AIR_DENSITY = 1.225;   ///< Standard air density at sea level (kg/m^3)
    
    /**
     * @brief Drag factor constant.
     * 
     * Calculated as: F = (rho * A) / (2 * m)
     * Area = 0.00143 m^2, F = 0.01907
     */
    constexpr double FACTOR_F = 0.01907; 
    
    // Drag Coefficients
    constexpr double CD_LAMINAR = 0.55;     ///< Drag coefficient for laminar flow (low speed)
    constexpr double CD_TURBULENT = 0.275;  ///< Drag coefficient for turbulent flow (high speed)
    
    // Transition Velocities (m/s)
    constexpr double V_TRANSITION_START = 18.0; ///< Velocity where laminar-turbulent transition begins
    constexpr double V_TRANSITION_END = 25.0;   ///< Velocity where transition completes (fully turbulent)
}

#endif // CONSTANTS_HPP