#ifndef TYPES_HPP
#define TYPES_HPP

#include <array>

// State vector indices
enum StateIndex { X_POS = 0, Y_POS = 1, X_VEL = 2, Y_VEL = 3 };

// 4D state: [x_pos, y_pos, x_vel, y_vel]
using State = std::array<double, 4>;

#endif // TYPES_HPP