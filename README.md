# ProjectileTrajectory

Compute the optimal launch angle for a projectile to achieve the furthest horizontal distance before striking the ground.

## Description


There exists an optimum angle $\theta_{opt}$ when launching a projectile subject to gravity, such that a maximum horizontal distance is flown before falling back to a reference launch height (or ground level).   This angle can be predicted analytically in the absence of atmospheric drag using the following formula.

$$\theta_{opt} = \arctan\left(\frac{v_0}{\sqrt{2gh + v_0^2}}\right)$$

Where $v_0$ represents initial velocity, $g$ is the gravitational acceleration constant, and $h$ is the launch height.  Notice a familiar result from elementary physics, with $h = 0$, $\theta_{opt} = 45^\circ$

In the presence of atmospheric drag, the optimal launch angle decreases. We don't have an analytic formula to calculate this result. Therefore we will simulate the atmospheric drag effect on the projectile.  Where drag is proportional to instantaneous velocity squared in each time step according to the following formula.

$$F_d = \frac{1}{2} \rho v^2 C_d A$$

$F_d = ma$ represents the drag force acting on the projectile. Consolidating factors that are assumed constant during the simulation, we'll parameterize the effect of drag and mass into a derivative function where the only tunable parameter is the inverse ballistic coefficient $\frac{k}{m}$ or "k_over_m".

[ballistic coefficient (m/Cd/A)](https://en.wikipedia.org/wiki/Ballistic_coefficient)

## Implementation Details

### Physics Model
The simulation solves the 2D projectile motion equations with quadratic air drag.
- **State Vector**: $\mathbf{S} = [x, y, v_x, v_y]$
- **Derivatives**:
  - $\dot{x} = v_x$
  - $\dot{y} = v_y$
  - $\dot{v_x} = -\frac{k}{m} v v_x$
  - $\dot{v_y} = -g - \frac{k}{m} v v_y$
  where $v = \sqrt{v_x^2 + v_y^2}$ and $\frac{k}{m}$ is the drag factor derived from $\frac{1}{2m} \rho C_d A$.

### Numerical Methods
- **Integration**: Fourth-order Runge-Kutta (RK4) is used to advance the state vector in discrete time steps ($\Delta t$). This offers a good balance between speed and accuracy compared to Euler or Heun methods.
- **Event Detection**: The simulation loop terminates when the vertical position $y < 0$. To determine the exact impact location, the program performs a linear interpolation between the last state above ground and the current state below ground.

### Optimization Strategy
To find the launch angle $\theta$ that maximizes horizontal distance $x_{impact}$, the program treats the simulation as a function $f(\theta) \rightarrow x_{impact}$.
- **Algorithm**: Golden Section Search.
- **Why**: This is a robust bracket-reduction method that finds extrema without requiring derivatives of the objective function. It iteratively narrows down the range $[a, b]$ containing the optimal angle until the width is within a specified tolerance.

## Build and Run

### Requirements
- C++17 compliant compiler
- CMake 3.10+

### Step-by-Step
1. **Clone the repository**:
   ```sh
   git clone <repository_url>
   cd ProjectileTrajectory
   ```

2. **Configure and Build**:
   ```sh
   mkdir build
   cd build
   cmake ..
   make
   ```

3. **Run**:
   ```sh
   ./main
   ```



## Getting Started

### Dependencies

* c++ standard version 14 compiler

### Building

A solution project can be built and run using the included CMakeLists.txt by executing the following commands (after installing cmake).  
```
mkdir build && cd build
cmake ..
cmake --build .
```

### Executing program

Run the executable `./main`

### Sample Result

```
Using drag coefficient k/m = 0.0025
Target distance precision: 0.1 m

Phase 1: Finding approximate optimal angle (coarse search):
Approximate optimal angle: 37.907528 degrees
Maximum distance: 398.029005 m

Phase 2: Finding angle precision for distance tolerance 0.100000 m:

Phase 3: Refining optimal angle to required precision:

=== RESULTS ===
Optimal launch angle: 37.912268 degrees
Maximum distance: 398.029002 m
Required angle precision: ±0.8 degrees
  (to maintain distance within ±0.100000 m)

Verification at angle tolerance boundaries:
  Angle: 37.149193° → Distance: 397.930246 m (Δ = 0.098756 m)
  Angle: 37.912268° → Distance: 398.029002 m (optimal)
  Angle: 38.675343° → Distance: 397.928587 m (Δ = 0.100415 m)
```

## Help

Should cmake not be available on the build machine. Compiling by command line on a mac is performed as follows.

`clang++ -std=c++14 -O3 -o main src/*.cpp -I include`

## Authors

Contributors names and contact info

Troy Wigton  

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MIT License - see the LICENSE.md file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [Ballistic_Coefficient](https://en.wikipedia.org/wiki/Ballistic_coefficient)
* [Markdown Preview](https://markdownlivepreview.com/)
* [Runge-Kutta Methods](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
