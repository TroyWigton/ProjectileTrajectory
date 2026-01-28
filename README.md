# ProjectileTrajectory

Compute the optimal launch angle for a projectile to achieve the furthest horizontal distance before striking the ground.

## Description


There exists an optimum angle $\theta_{opt}$ when launching a projectile subject to gravity, such that a maximum horizontal distance is flown before falling back to a reference launch height (or ground level).   This angle can be predicted analytically in the absence of atmospheric drag using the following formula.

$$\theta_{opt} = \arctan\left(\frac{v_0}{\sqrt{2gh + v_0^2}}\right)$$

Where $v_0$ represents initial velocity, $g$ is the gravitational acceleration constant, and $h$ is the launch height.  Notice a familiar result from elementary physics, with $h = 0$, $\theta_{opt} = 45^\circ$

In the presence of atmospheric drag, the optimal launch angle decreases. We don't have an analytic formula to calculate this result. Therefore we will simulate the atmospheric drag effect on the projectile.  Where drag is proportional to instantaneous velocity squared in each time step according to the following formula.

$$F_d = \frac{1}{2} \rho v^2 C_d A$$

$F_d = ma$ represents the drag force acting on the projectile. Consolidating factors that are assumed constant during the simulation, we'll parameterize the effect of drag and mass within the derivative function as `k_over_m`,  $\frac{k}{m} = \frac{\rho C_d A}{2m}$.

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

#### Example Coefficients

Below are calculated `k_over_m` values for common objects, assuming standard sea-level air density ($\rho \approx 1.225 \text{ kg/m}^3$).

| Object | Mass (kg) | Diam (m) | Cd | Area (m^2) | k/m (1/m) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Vacuum** (No Drag) | - | - | 0.0 | - | **0.0** |
| **9mm Bullet** (115gr FMJ) | 0.00745 | 0.009 | 0.295 | $6.36 \times 10^{-5}$ | **0.00154** |
| **Golf Ball** | 0.04593 | 0.0427 | 0.30 | $1.43 \times 10^{-3}$ | **0.00573** |
| **Ping Pong Ball** | 0.0027 | 0.040 | 0.47 | $1.26 \times 10^{-3}$ | **0.13398** |

*Note: These are approximations using constant drag coefficients. Real-world ballistics involve velocity-dependent $C_d$ (Mach number variations), but these values provide good baselines for comparison.*

### Numerical Methods
- **Integration**: Fourth-order Runge-Kutta (RK4) is used to advance the state vector in discrete time steps ($\Delta t$). This offers a good balance between speed and accuracy compared to Euler or Heun methods.
- **Event Detection**: The simulation loop terminates when the vertical position $y < 0$. To determine the exact impact location, the program performs a linear interpolation between the last state above ground and the current state below ground.

### Optimization Strategy
To find the launch angle $\theta$ that maximizes horizontal distance $x_{impact}$, the program treats the simulation as a function $f(\theta) \rightarrow x_{impact}$.
- **Algorithm**: Golden Section Search.
- **Why**: This is a robust bracket-reduction method that finds extrema without requiring derivatives of the objective function. It iteratively narrows down the range $[a, b]$ containing the optimal angle until the width is within a specified tolerance.

## Build and Run

### Requirements
- C++14 compliant compiler
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

3. **Run Main Simulation**:
   ```sh
   ./main
   ```

   **Sample Result**:
   ```
   Using drag coefficient k/m = 0.0057
   Target distance precision: 0.01 m

   Phase 1: Finding approximate optimal angle (coarse search):
   Approximate optimal angle: 34.735268 degrees
   Maximum distance: 246.274544 m

   Phase 2: Finding angle precision for distance tolerance 0.010000 m:

   Phase 3: Refining optimal angle to required precision:

   === RESULTS ===
   Optimal launch angle: 34.737258 degrees
   Maximum distance: 246.274544 m
   Required angle precision: ±0.3 degrees
     (to maintain distance within ±0.010000 m)

   Verification at angle tolerance boundaries:
     Angle: 34.416829° → Distance: 246.264668 m (Δ = 0.009876 m)
     Angle: 34.737258° → Distance: 246.274544 m (optimal)
     Angle: 35.057687° → Distance: 246.264531 m (Δ = 0.010013 m)
   ```

4. **Run Integrator Comparison**:
   To benchmark different numerical methods (Euler, Heun, RK4, RK8) against a high-precision ground truth:
   ```
   ./compare_integrators
   ```
   This tool iteratively refines the time step ($\Delta t$) for each integrator until it matches the ground truth distance within a specified tolerance (0.01m), reporting the efficiency (steps required) of each method.

   **Sample Result**:
   ```
   Project Comparison Tool
   Drag coefficient k/m = 0.134
   Target precision: 0.0001 m

   Optimal angle used for comparison: 24.3629 degrees
   Ground Truth (RK8, dt=1e-05): 22.85398 m

   Comparison: Required Steps & dt to reach 0.00010 m precision
   --------------------------------------------------------------------------------
   Integrator     Steps          dt (s)         Distance (m)        Error (m)
   --------------------------------------------------------------------------------
   Euler          1073219        1.90735e-06    22.85390            8.25e-05
   Heun           4193           0.000488281    22.85401            2.57e-05
   RK4            263            0.0078125      22.85397            1.20e-05
   RK8            132            0.015625       22.85398            3.39e-06
   ```

5. **Run Unit Tests**:
   Executes a suite of validation tests to ensure physics engine accuracy.
   - **Vacuum Test**: Verifies $\theta_{opt} \approx 45^\circ$ for all integrators when drag is zero.
   - **Consistency Test**: Cross-references lower-order integrators (Euler, Heun) against the high-precision RK8 method to ensure results are within expected convergence limits under drag conditions.
   ```sh
   ./test
   ```
   **Sample Results**
   ```
   Starting Test Suite: Vacuum Trajectory Optimization
   Target Angle: 45 degrees
   Tolerance: 0.05 degrees
   k/m: 0
   
   Testing Integrator: Euler, Derivative: V Squared Drag... PASSED (Optimal Angle: 44.9998)
   Testing Integrator: Euler, Derivative: Linear Drag... PASSED (Optimal Angle: 44.9998)
   Testing Integrator: Euler, Derivative: No Drag... PASSED (Optimal Angle: 44.9998)
   Testing Integrator: Heun, Derivative: V Squared Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: Heun, Derivative: Linear Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: Heun, Derivative: No Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: RK4, Derivative: V Squared Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: RK4, Derivative: Linear Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: RK4, Derivative: No Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: RK8, Derivative: V Squared Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: RK8, Derivative: Linear Drag... PASSED (Optimal Angle: 44.9999)
   Testing Integrator: RK8, Derivative: No Drag... PASSED (Optimal Angle: 44.9999)
   
   --------------------------------------------------
   Starting Test Suite: Integrator Consistency (Non-zero Drag)
   k/m: 0.0057
   
   Evaluating Derivative Function: V Squared Drag
          Euler: Angle = 34.7324, Distance = 246.264
           Heun: Angle = 34.7358, Distance = 246.275
            RK4: Angle = 34.7358, Distance = 246.275
            RK8: Angle = 34.7358, Distance = 246.275
     -> Consistency Check (Euler vs RK8): PASSED
        (Diffs - Angle: 0.00343658, Dist: 0.0109617)
     -> Consistency Check (Heun vs RK8): PASSED
        (Diffs - Angle: 0, Dist: 3.66264e-08)
     -> Consistency Check (RK4 vs RK8): PASSED
        (Diffs - Angle: 0, Dist: 9.37916e-13)
   
   Evaluating Derivative Function: Linear Drag
          Euler: Angle = 44.2345, Distance = 966.492
           Heun: Angle = 44.2359, Distance = 966.425
            RK4: Angle = 44.2359, Distance = 966.425
            RK8: Angle = 44.2359, Distance = 966.425
     -> Consistency Check (Euler vs RK8): PASSED
        (Diffs - Angle: 0.00141116, Dist: 0.0670337)
     -> Consistency Check (Heun vs RK8): PASSED
        (Diffs - Angle: 0, Dist: 1.27375e-07)
     -> Consistency Check (RK4 vs RK8): PASSED
        (Diffs - Angle: 0, Dist: 1.80762e-11)
   
   All tests passed successfully!
   ```

6. **Run Drag Coefficient Analysis**:
   Demonstrates how the optimal launch angle and maximum distance vary with increasing drag ($\frac{k}{m}$ from 0.0 to 0.2).
   ```sh
   ./compare_k_over_m
   ```
   **Sample Results**:
   ```
   K/M Variation Analysis
   --------------------------------------------------------
   v0 = 100 m/s
   integrator = RK4, dt = 0.001 s
   --------------------------------------------------------
          k/m      Optimal Angle (deg)    Max Distance (m)
   --------------------------------------------------------
      0.00000                 44.99986          1019.36799
      0.00100                 40.92337           600.84063
      0.00124                 40.30281           553.28964
      0.00153                 39.63102           505.50411
      0.00189                 38.91255           458.30923
      0.00233                 38.16013           412.46152
      0.00289                 37.37236           368.60698
      0.00357                 36.56577           327.25561
      0.00441                 35.74291           288.77216
      0.00545                 34.91226           253.38064
      0.00674                 34.08320           221.17880
      0.00833                 33.25647           192.15838
      0.01029                 32.44274           166.22760
      0.01272                 31.64991           143.23317
      0.01572                 30.86820           122.98019
      0.01943                 30.11288           105.24900
      0.02402                 29.38003            89.80873
      0.02969                 28.66842            76.42763
      0.03670                 27.98629            64.88067
      0.04537                 27.33439            54.95465
      0.05608                 26.69778            46.45142
      0.06931                 26.09051            39.18970
      0.08568                 25.50826            33.00564
      0.10590                 24.95898            27.75276
      0.13090                 24.42804            23.30121
      0.16180                 23.91749            19.53672
      0.20000                 23.42878            16.35945
      ```

## Help

If CMake is not available, you can compile the individual programs manually using a C++14 compiler.

**Main Simulation:**
```sh
clang++ -std=c++14 -O3 -I include -o main src/main.cpp src/simulation.cpp src/derivative_functions.cpp
```

**Integrator Comparison:**
```sh
clang++ -std=c++14 -O3 -I include -o compare_integrators src/compare_integrators.cpp src/simulation.cpp src/derivative_functions.cpp
```

**Unit Tests:**
```sh
clang++ -std=c++14 -O3 -I include -o test src/test.cpp src/simulation.cpp src/derivative_functions.cpp
```

**Drag Coefficient Analysis:**
```sh
clang++ -std=c++14 -O3 -I include -o compare_k_over_m src/compare_k_over_m.cpp src/simulation.cpp src/derivative_functions.cpp
```

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
