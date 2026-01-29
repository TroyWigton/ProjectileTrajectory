#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

// Physical constants
constexpr double GRAVITY_EARTH = 9.81;  // m/s²

// Drag coefficients (k/m) references
namespace DragRatios {
    constexpr double VACUUM = 0.0;
    constexpr double BULLET_223 = 0.0012; // Approx for .223 Rem (55gr), derived from Cd ~ 0.23
    constexpr double BULLET_9MM = 0.0015;
    constexpr double GOLF_BALL = 0.0057;
    constexpr double PING_PONG_BALL = 0.134;
}

namespace GolfBallPhysics {
    // Physical properties
    constexpr double MASS_KG = 0.04593;
    constexpr double DIAMETER_M = 0.04267;
    constexpr double AIR_DENSITY = 1.225;
    
    // Constant Factor F = (rho * A) / (2 * m)
    // Area = 0.00143 m^2, F = 0.01907
    constexpr double FACTOR_F = 0.01907; 
    
    // Drag Coefficients
    constexpr double CD_LAMINAR = 0.55;
    constexpr double CD_TURBULENT = 0.275;
    
    // Transition Velocities (m/s)
    constexpr double V_TRANSITION_START = 18.0;
    constexpr double V_TRANSITION_END = 25.0;
}

#endif // CONSTANTS_HPP