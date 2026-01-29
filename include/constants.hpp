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

#endif // CONSTANTS_HPP