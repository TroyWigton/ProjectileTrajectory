#include<array>
#include<functional>

// 4D state: [x_pos, y_pos, x_vel, y_vel]
typedef std::array<double,4> State;

template<typename State>
using IntegratorFunc = State(*)(const State&, double, double, 
                               std::function<void(const State&, double, State&)>);

//Forward Euler (Explicit Euler)
template<typename DerivativeFunc>
State euler_step(const State& state, double t, double dt, DerivativeFunc deriv_func);

//Heun's Method (Trapezoidal Rule)
template<typename DerivativeFunc>
State heun_step(const State& state, double t, double dt, DerivativeFunc deriv_func);

//Fourth Order Runge-Kutta
template<typename DerivativeFunc>
State rk4_step(const State& state, double t, double dt, DerivativeFunc deriv_func);

// 8th Order Runge-Kutta
template<typename State>
State rk8_step(const State& state, double t, double dt, std::function<void(const State&, double, State&)> deriv_func);
