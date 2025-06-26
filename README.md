# ProjectileTrajectory

Compute the optimal launch angle for a projectile to achieve the furthest horizontal distance before striking the ground.

## Description

We learned in physics class that launching a projectile at 45 degrees from horizontal achieves the optimal horizontal distance travled before the object returns to the launch height level (ground).  We learned how to prove this result in the absense of any aerodynamic drag.  But questions remain. Is 45 degress still an optimal angle when launched from above ground level?  Does the optimum solution differ in the presense of aerodynamic drag, and if so, by what measure, based on representative factors such as [balistic coefficient (m/Cd/A)](https://en.wikipedia.org/wiki/Ballistic_coefficient) or its inverse k/m.

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

`clang++ -std=c++14 -o3 -o main src/*.cpp -I/include`

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
