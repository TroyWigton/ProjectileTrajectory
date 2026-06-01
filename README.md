# ProjectileTrajectory

Compute the optimal launch angle for a projectile to achieve the furthest horizontal distance before striking the ground.

## Contents

- [The Problem](#the-problem)
- [Numerical Methods](#numerical-methods)
  - [Integrators](#integrators)
  - [Golden Section Search](#golden-section-search)
  - [Event Detection](#event-detection)
- [Physics Models](#physics-models)
- [Getting Started](#getting-started)
- [Executables](#executables)
  - [optimal_angle_search](#optimal_angle_search)
  - [compare_integrators](#compare_integrators)
  - [convergence_study](#convergence_study)
  - [compare_k_over_m](#compare_k_over_m)
  - [mach_drag_study](#mach_drag_study)
  - [golf_ball_drag_study](#golf_ball_drag_study)
  - [demo_rk45_adaptive](#demo_rk45_adaptive)
  - [example_3d_motion](#example_3d_motion)
- [Validation and Testing](#validation-and-testing)
- [Manual Compilation](#manual-compilation)
- [Authors](#authors)
- [Version History](#version-history)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## The Problem

There exists an optimum angle $\theta_{opt}$ when launching a projectile subject to gravity, such that a maximum horizontal distance is flown before falling back to a reference launch height (or ground level).   This angle can be predicted analytically in the absence of atmospheric drag using the following formula.

$$\theta_{opt} = \arctan\left(\frac{v_0}{\sqrt{2gh + v_0^2}}\right)$$

The corresponding maximum horizontal distance $R_{max}$ achieved at this angle is:

$$R_{max} = \frac{v_0}{g}\sqrt{2gh + v_0^2}$$

Where $v_0$ represents initial velocity, $g$ is the gravitational acceleration constant, and $h$ is the launch height.  Notice a familiar result from elementary physics, with $h = 0$, $\theta_{opt} = 45^\circ$

In the presence of atmospheric drag, the optimal launch angle decreases. We don't have an analytic formula to calculate this result. Therefore we will simulate the atmospheric drag effect on the projectile.  Where drag is proportional to instantaneous velocity squared in each time step according to the following formula.

$$F_d = \frac{1}{2} \rho v^2 C_d A$$

$F_d = ma$ represents the drag force acting on the projectile. Consolidating factors that are assumed constant during the simulation, we'll parameterize the effect of drag and mass within the derivative function as `k_over_m`,  $\frac{k}{m} = \frac{\rho C_d A}{2m}$.

## Numerical Methods

### Integrators

The integrators live in `include/integrators.hpp` as **header-only template functions, generic over both the state type and the derivative function**. Each method shares the uniform signature `(state, t, dt, deriv_func)` and returns the next state, so the same code integrates the 4-DOF projectile state `[x, y, vx, vy]`, the 6-DOF `State6D`, or any other state type without modification. Stage temporaries are `thread_local`, which keeps `std::array`-backed states allocation-free across millions of steps. Six methods are provided:

- **Euler (1st Order)**: Basic forward stepping method used as a baseline. Simple but requires extremely small time steps for accuracy.
- **Heun (2nd Order)**: Also known as the trapezoidal method. Offers improved stability over Euler.
- **RK4 (4th Order)**: The classic Runge-Kutta method. Provides an excellent balance of speed and accuracy; used as the default engine.
- **RK45 (Dormand-Prince)**: An embedded 5(4) method available in two variants — `rk45_step` (fixed dt, high-precision) and `rk45_adaptive_step` (returns an error estimate alongside the new state for step-size control). Both share the same 7-stage tableau.
- **RK8 (8th Order)**: High-order method requiring 13 function evaluations per step. Extremely accurate, used as the "Ground Truth" for benchmarking other methods.

### Golden Section Search

To find the launch angle $\theta$ that maximizes horizontal distance $x_{impact}$, the program treats the simulation as a function $f(\theta) \rightarrow x_{impact}$. The search itself is a template, `golden_section_search_max(func, a, b, tol)` in `include/golden_section.hpp`, generic over the objective callable — `optimal_angle_search` simply passes it a lambda that wraps `Simulation::run`.

- **Algorithm**: Golden Section Search.
- **Why**: This is a robust bracket-reduction method that finds extrema without requiring derivatives of the objective function. It iteratively narrows down the range $[a, b]$ containing the optimal angle until the width is within a specified tolerance.
- **Process**:
    1.  **Initialization**: Start with a known bracket $[a,b]$ (e.g., $10^\circ$ to $80^\circ$) that is certain to contain the maxima.
    2.  **Evaluation**: Evaluate the distance function at two interior points $c$ and $d$ chosen by the golden ratio $\phi = \frac{1+\sqrt{5}}{2}$.
    3.  **Reduction**: Discard the sub-segment that does not contain the maximum, maintaining the golden ratio property for the next iteration.
    4.  **Termination**: Repeat until $(b-a)$ is less than the desired angular precision.
- **Refinement**: To improve efficiency, the program employs a multi-phase approach. It first performs a coarse search to locate the approximate peak, calculates the sensitivity of distance to angle near that peak, and then performs a high-precision refinement focused narrowly on the optimal region.

### Event Detection

The main simulation loop terminates when the vertical position $y < 0$. To determine the exact impact location, the program performs a linear interpolation between the last state above ground and the current state below ground.

## Physics Models

The simulation solves the 2D projectile motion equations with quadratic air drag.
- **State Vector**: $\mathbf{S} = [x, y, v_x, v_y]$
- **Derivatives**:
  - $\dot{x} = v_x$
  - $\dot{y} = v_y$
  - $\dot{v_x} = -\frac{k}{m} v v_x$
  - $\dot{v_y} = -g - \frac{k}{m} v v_y$ , where $v = \sqrt{v_x^2 + v_y^2}$ and $\frac{k}{m} = \frac{\rho C_d A}{2m}$

### Example Coefficients

Below are calculated `k_over_m` values for common objects, assuming standard sea-level air density ($\rho \approx 1.225 \text{ kg/m}^3$).

| Object | Mass (kg) | Diam (m) | Cd | Area (m^2) | k/m (1/m) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Vacuum** (No Drag) | - | - | 0.0 | - | **0.0** |
| **9mm Bullet** (115gr FMJ) | 0.00745 | 0.009 | 0.295 | $6.36 \times 10^{-5}$ | **0.00154** |
| **Golf Ball** | 0.04593 | 0.0427 | 0.30 | $1.43 \times 10^{-3}$ | **0.00573** |
| **Ping Pong Ball** | 0.0027 | 0.040 | 0.47 | $1.26 \times 10^{-3}$ | **0.13398** |

**Note:** These are approximations using constant drag coefficients. Real-world ballistics involve velocity-dependent $C_d$ (Mach number and Reynolds number variations), but these values provide sufficient baselines for comparison.

## Getting Started

### Requirements
- C++14 compliant compiler
- CMake 3.10+

### Build

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
   cmake --build . -j
   ```

All targets are built to `build/`. Run any of them from there (e.g. `./optimal_angle_search`); see [Executables](#executables) for what each one does. If CMake is unavailable, see [Manual Compilation](#manual-compilation).

## Executables

### optimal_angle_search

Finds the optimal launch angle for maximum distance via golden section search.

```sh
./optimal_angle_search
```

<details>
<summary>Sample output</summary>

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

</details>

### compare_integrators

Benchmarks different numerical methods (Euler, Heun, RK4, RK45, RK8) by measuring the computational effort required to achieve high accuracy.
- **Scenario**: A Ping Pong ball (High Drag, $k/m \approx 0.134$) launched at $200$ m/s at a fixed $45^\circ$ angle.
- **Method**: The tool runs a fixed-time simulation ($T=20.0s$) and compares the final position against a "Ground Truth" calculated using RK8 with an extremely fine time step ($\Delta t = 10^{-5}$ s).
- **Optimization**: It iteratively adjusts the time step ($\Delta t$) for each integrator until its final position error is within the target precision ($0.001$ m).
- **Output**: Reports the number of steps, time step, final position error, and wall-clock runtime (one timed run at the converged dt) for each method, illustrating the efficiency trade-off between lower-order methods (millions of cheap steps) and higher-order methods (fewer but more expensive steps).

```sh
./compare_integrators
```

<details>
<summary>Sample output</summary>

```
Numerical Integrator Benchmark
Simulates a 2D projectile (gravity + v^2 air drag) for a fixed duration T,
then compares each integrator's final (x, y) position against an RK8 ground
truth to measure accuracy vs. computational cost.

Launch Angle:         45 deg
Initial Velocity (v0): 200 m/s
Drag/Mass Ratio (k/m): 0.134
Simulation Duration:  20 s
Target Precision:     0.001 m (Position Error Norm at t=T)

Reference Truth (RK8, dt=1e-05s) at t=20s: (x=24.564 m, y=-141.258 m)

Methodology: Iteratively reducing time step (dt) until position error < target precision.
             Runtime is measured by re-running once at the converged dt.
Performance Comparison: Steps, dt, error, and wall-clock runtime to meet target precision
--------------------------------------------------------------------------------------------
Integrator     Steps          dt (s)         Final Error (m)     Runtime (ms)
--------------------------------------------------------------------------------------------
Euler          2156009        9.2764e-06     9.99e-04            50.722
Heun           14309          0.00139772     9.88e-04            0.682
RK4            1787           0.0111919      9.88e-04            0.157
RK45           1568           0.0127551      9.48e-04            0.205
RK8            429            0.04662        9.83e-04            0.120
```

</details>

### convergence_study

Measures the empirical order of convergence for each fixed-step integrator.
For each method, the same smooth v² drag projectile is integrated to a fixed
t=T at progressively halved dt values; the global position error at t=T is
compared against an ultra-precise RK8 reference, and the empirical order
$= \log_2(\text{error}_{prev} / \text{error}_{curr})$ is reported per row.
For a method of theoretical order $p$, the empirical order should approach
$p$ until the error reaches the floating-point precision floor.

```sh
./convergence_study
```

<details>
<summary>Sample output</summary>

```
Numerical Integrator Convergence Test
------------------------------------------------------------
Test problem (smooth v^2 drag):
  v0           = 100 m/s
  launch angle = 45 deg
  k/m          = 0.0057 (golf ball)
  duration T   = 5 s
Reference: RK8 at dt = 1e-06 (truth = (176.396, 88.2665) m)

Each integrator is run at progressively halved dt values; the
global position error at t=T is computed against the reference.
Empirical order = log2(error_prev / error_curr). For an order-p
method, the empirical order should approach p until the error
reaches the floating-point precision floor (~1e-13).

------------------------------------------------------------
Euler (theoretical order 1)
------------------------------------------------------------
dt (s)      Error (m)      Error ratio    Emp. order
0.1         4.507e-01      -              -
0.05        2.220e-01      2.03           1.02
0.025       1.104e-01      2.01           1.01
0.0125      5.506e-02      2              1.00
0.00625     2.750e-02      2              1.00
0.003125    1.374e-02      2              1.00

------------------------------------------------------------
Heun (theoretical order 2)
------------------------------------------------------------
dt (s)      Error (m)      Error ratio    Emp. order
0.1         1.953e-02      -              -
0.05        4.495e-03      4.34           2.12
0.025       1.079e-03      4.16           2.06
0.0125      2.646e-04      4.08           2.03
0.00625     6.549e-05      4.04           2.01
0.003125    1.629e-05      4.02           2.01

------------------------------------------------------------
RK4 (theoretical order 4)
------------------------------------------------------------
dt (s)      Error (m)      Error ratio    Emp. order
0.1         3.088e-05      -              -
0.05        1.828e-06      16.9           4.08
0.025       1.112e-07      16.4           4.04
0.0125      6.871e-09      16.2           4.02
0.00625     4.483e-10      15.3           3.94

------------------------------------------------------------
RK45 (fixed) (theoretical order 5)
------------------------------------------------------------
dt (s)      Error (m)      Error ratio    Emp. order
0.5         5.390e-03      -              -
0.25        4.752e-05      113            6.83
0.125       4.723e-07      101            6.65
0.0625      4.209e-09      112            6.81
0.03125     7.556e-11      55.7           5.80

------------------------------------------------------------
RK8 (theoretical order 8)
------------------------------------------------------------
dt (s)      Error (m)      Error ratio    Emp. order
1           2.136e-05      -              -
0.5         7.318e-08      292            8.19
0.25        2.984e-10      245            7.94
0.125       8.228e-11      3.63           1.86
0.0625      8.229e-11      1              -0.00
```

</details>

Euler, Heun, and RK4 each recover their theoretical orders cleanly across
the dt sweep. RK45 shows super-convergence at moderate dt (next-order error
terms partially cancel the leading term), with the empirical order trending
toward the theoretical 5 as dt shrinks. RK8 recovers order 8 cleanly for
two halvings before the error reaches the floating-point precision floor
(~1e-10 for this problem, limited by the reference's own accumulated
rounding).

### compare_k_over_m

Demonstrates how the optimal launch angle and maximum distance vary with increasing drag ($\frac{k}{m}$ from 0.0 to 0.2).

```sh
./compare_k_over_m
```

<details>
<summary>Sample output</summary>

```
K/M Variation Analysis
--------------------------------------------------------
Sweeps the drag coefficient (k/m) over 0 to 0.2 and finds the launch angle
that maximizes ground-impact distance for each k/m (golden section search).
Each trial integrates the 2D trajectory until y < 0, then interpolates to
the exact ground crossing.

Drag Model:            v^2 quadratic drag (drag_deriv_v_squared)
                       a_drag = -(k/m) * |v| * v_vec
Integrator:            RK4, fixed dt = 0.001 s
Initial Velocity (v0): 100 m/s
Launch Height (h0):    0 m
Gravity (g):           9.81 m/s^2
Angle Search Bracket:  [10, 80] deg, tolerance 1e-4 deg
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

</details>

### mach_drag_study

Evaluates the impact of Mach-dependent drag (transonic wave drag with supersonic fade) on a high-velocity projectile (.223 Remington, ~Mach 2.84 at the muzzle).
It compares a vacuum baseline, a v^2 drag model with constant $C_d$, and a v^2 drag model with variable $C_d$ that rises through the transonic regime and fades back at higher Mach.

```sh
./mach_drag_study
```

<details>
<summary>Sample output</summary>

```
Supersonic Variable Drag Model Evaluation (.223 Remington)
--------------------------------------------------------
Initial Velocity: 975 m/s (Mach 2.8426)
Base k/m: 0.0012
Time step: 0.0010 s
--------------------------------------------------------

Model 0: Vacuum (Theoretical Baseline)
  Optimal Angle: 44.9999 degrees
  Max Distance:  96903.67 m

Model 1: v^2 drag, constant Cd
  Optimal Angle: 24.7662 degrees
  Max Distance:  2481.73 m

Model 2: v^2 drag, variable Cd (transonic rise with supersonic fade)
  Optimal Angle: 26.3674 degrees
  Max Distance:  2162.62 m

Impact Analysis (Model 2 variable Cd vs Model 1 constant Cd):
  Range Reduction: 319.10 m (12.86% loss vs constant Cd)
  Angle Shift:     1.60 degrees (variable Cd - constant Cd)
```

</details>

### golf_ball_drag_study

Compares optimal-angle trajectories under three golf-ball drag models: vacuum, constant-Cd v² drag, and velocity-dependent Cd that transitions from laminar to turbulent across the 18-25 m/s range.

```sh
./golf_ball_drag_study
```

<details>
<summary>Sample output</summary>

```
Golf Ball Trajectory Optimization (v0=70 m/s ~ 156 mph ball speed)
------------------------------------------------------------------------
Cd model parameters:
  Constant Cd: 0.4125  (mean of laminar=0.550 and turbulent=0.275)
  Variable Cd: 0.550 (laminar, v < 18.0 m/s) -> 0.275 (turbulent, v >= 25.0 m/s)
               with linear transition between 18.0 and 25.0 m/s
------------------------------------------------------------------------
1. Vacuum                  Optimal Angle: 44.9987 deg | Max Distance: 499.490 m
2. v^2 drag, constant Cd   Optimal Angle: 36.2675 deg | Max Distance: 153.276 m
3. v^2 drag, variable Cd   Optimal Angle: 37.8064 deg | Max Distance: 192.229 m
```

</details>

### demo_rk45_adaptive

Demonstrates the adaptive RK45 integrator on a golf-ball trajectory at two error tolerances. Reports accepted/rejected step counts, achieved min/max dt range, and final flight distance.

```sh
./demo_rk45_adaptive
```

### example_3d_motion

Standalone 3D-projectile demo using the 6-DOF `State6D = [x, y, z, vx, vy, vz]`. Confirms that the integrator templates accept arbitrary state dimensions, not just the 4-DOF state used by the projectile motion code.

```sh
./example_3d_motion
```

<details>
<summary>Sample output</summary>

```
Running 3D Projectile Simulation (6-DOF State)...
Drag (k/m): 0.005
t=1.0000 x=56.7092 y=53.2586 z=9.9994 v=61.9955
t=2.0000 x=97.8509 y=83.3161 z=17.2538 v=42.8891
Impact at t=9.0200 x=253.9047 y=-0.2314 z=44.7703
```

</details>

## Validation and Testing

The `test` executable runs a suite of validation tests to ensure physics-engine accuracy.
- **Vacuum Test**: Verifies $\theta_{opt} \approx 45^\circ$ for all integrators when drag is zero.
- **Consistency Test**: Cross-references lower-order integrators (Euler, Heun, RK4) against the high-precision RK8 method to ensure results are within expected convergence limits under drag conditions.

```sh
./test
```

<details>
<summary>Sample output</summary>

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

--------------------------------------------------
Starting Test Suite: Vacuum Trajectory with Height (No Drag)
Launch Height h0: 100 m
Initial Velocity v0: 100 m/s
Gravity g: 9.81 m/s^2
Theoretical Optimal Angle: 42.4373 degrees
Experimental Optimal Angle: 42.4377 degrees
PASSED (Difference: 0.000363027)
Theoretical Max Distance: 1114.89 m
Experimental Max Distance: 1114.89 m
Distance Check: PASSED (Difference: 1.02653e-06)

--------------------------------------------------
Starting Test Suite: V Squared Drag Trajectory with Variable Height
Parameters: v0 = 100 m/s, k/m = 0.0057 (Golf Ball)
Hypothesis: As launch height increases, optimal angle should decrease.

Height (m) Optimal Angle (deg)    Max Distance (m)
--------------------------------------------------
         0             34.7358             246.275
        50             29.3616             269.451
       100              25.285             286.965
       200             19.7593             310.601
       500             13.4513              337.36
Trend verification: PASSED

All tests passed successfully!
```

</details>

## Manual Compilation

If CMake is not available, you can compile the individual programs manually using a C++14 compiler.

<details>
<summary>Per-target <code>clang++</code> commands</summary>

**Optimal-Angle Search:**
```sh
clang++ -std=c++14 -O3 -I include -o optimal_angle_search src/optimal_angle_search.cpp src/simulation.cpp src/derivative_functions.cpp
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

**Supersonic Drag Study:**
```sh
clang++ -std=c++14 -O3 -I include -o mach_drag_study src/mach_drag_study.cpp src/simulation.cpp src/derivative_functions.cpp
```

**Golf Ball Drag Study:**
```sh
clang++ -std=c++14 -O3 -I include -o golf_ball_drag_study src/golf_ball_drag_study.cpp src/simulation.cpp src/derivative_functions.cpp
```

**Adaptive RK45 Demo:**
```sh
clang++ -std=c++14 -O3 -I include -o demo_rk45_adaptive src/demo_rk45_adaptive.cpp src/simulation.cpp src/derivative_functions.cpp
```

**Convergence Study:**
```sh
clang++ -std=c++14 -O3 -I include -o convergence_study src/convergence_study.cpp src/simulation.cpp src/derivative_functions.cpp
```

**3D Motion Demo:**
```sh
clang++ -std=c++14 -O3 -I include -o example_3d_motion src/example_3d_motion.cpp
```

</details>

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
* [Projectile Motion](https://en.wikipedia.org/wiki/Projectile_motion)
* [Golden-section search](https://en.wikipedia.org/wiki/Golden-section_search)
* [Markdown Preview](https://markdownlivepreview.com/)
* [Runge-Kutta Methods](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
