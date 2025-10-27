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
With k_over_m = 0.0025
Maximum distance without drag: 254.8 m at 45.0000 degrees
Maximum distance with drag: 175.3 m at 42.1080 degrees
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
