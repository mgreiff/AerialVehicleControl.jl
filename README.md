# AerialVehicleControl.jl
This software serves as
a framework for the analysis and generation of controllers pertaining to
quad-rotor dynamics. The controllers are implemented in C89, and utilize
[LAPACK](http://www.netlib.org/lapack/) and
[BLAS](http://www.netlib.org/blas/), or the
[CMSIS DSP ARM math library](http://www.keil.com/pack/doc/CMSIS/DSP/html/index.html),
should these be available - if not, the
code runs without any dependencies. The general idea is to call these
controllers through Julia, to (i) test all of the components of the controllers
using unit testing, (ii) quantitatively study their closed loop
behavior in simulation, and (iii) tune them and offline in a data-driven manner.

Importantly, once a controller has been chosen for a given problem, the exact
code used in the simulations can be run directly on the application, being
platform independent C-code. Thus the project does not aim to solve a single
control problem, but rather serves as an analysis tool and implementation aid
for any given aerial vehicle control problem.

### Installation
To run the code, you will need to install

* [Julia](https://julialang.org/downloads/platform/)
* [GCC](https://gcc.gnu.org/)
* [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/)
* [Doxygen](https://www.doxygen.nl/download.html)

The last point is only required if you wish to regenerate the docs (see the
above links for platform specific installation instructions).

Clone the repository and enter the Julia REPL by
```
git clone https://github.com/mgreiff/AerialVehicleControl.jl
cd AerialVehicleControl.jl
julia
```
and activate the package it by running the following command in the Julia REPL
```
]activate .
```
You should now be able to run the tests and examples /test and /examples
respectively, and the complete documentation can be found in an interactive
HTML-format here, as generated from docs/DOCS.md.

### Citation
Academic citations are welcomed, for the FSF and FOF controllers respectively:
```
@article{mgreiff2020robust,
  title={Attitude Control on SU(2): Stability, Robustness, and Similarities},
  author={Greiff, Marcus and Robertsson, Anders and Sun, Zhiyong},
  journal={IEEE control systems letters},
  year={2021}
}

@book{08166fa7ac4643149ea0b29bf5fe5d09,
title = "Filtered Output Feedback Tracking Control of a Quadrotor UAV",
author = "Erjen Lefeber and Marcus Greiff and Anders Robertsson",
year = "2020",
month = may,
day = "9",
language = "English",
series = "DC Reports",
publisher = "Eindhoven University of Technology, Dynamics and Control Group, Department of Mechanical Engineering, Eindhoven, The Netherlands",
number = "DC 2020.053",
}
```
The former paper is to appear in control systems letters in early 2021.
