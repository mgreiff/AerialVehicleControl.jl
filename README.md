# AerialVehicleControl.jl
This software serves as a framework for the generation, testing and analysis of
controllers pertaining to quad-rotor dynamics. The controllers are implemented
in C89, utlilizing [LAPACK](http://www.netlib.org/lapack/) and
[BLAS](http://www.netlib.org/blas/). Soon it will also be possible to run with the
[CMSIS DSP ARM math library](http://www.keil.com/pack/doc/CMSIS/DSP/html/index.html),
for implementation on MCUs should, and also standalone without any dependencies.

The general idea is to evaluate through a separate implementation in Julia, to
1. test all of the components of the controllers using unit testing,
2. quantitatively study their closed loop behavior in simulation,
3. tune them offline using semidefinite programming.

Importantly, once a controller has been chosen and tuned for a given control problem,
the exact code used in the simulations can be run directly on the application, being
platform independent C-code. Thus the project does not aim to solve a single
control problem, but rather serves as an analysis tool and implementation aid
for any given aerial vehicle control problem.

For a holistic overview of the codebase, see the HTML documentation
[here](http://www.control.lth.se/personnel/marcus-greiff/).

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
Note that the former reference is currently incomplete, as the paper is to
appear in control systems letters in early 2021.

### Installation and Usage
To run the code, you will need to install

* [Julia](https://julialang.org/downloads/platform/) (last run on 1.5.1, any >=1.4 should work)
* [GCC](https://gcc.gnu.org/) (last run with version 9.3.0)
* [LAPACK](http://www.netlib.org/lapack/) and [BLAS](http://www.netlib.org/blas/) (last run on version 3.9.0-1build1)
* [Doxygen](https://www.doxygen.nl/download.html) (last run on version 1.8.17)

The last point is only required if you wish to regenerate the docs (see the
above links for platform specific installation instructions).

Once installed, clone the repository and enter the Julia REPL by
```
git clone https://github.com/mgreiff/AerialVehicleControl.jl
cd AerialVehicleControl.jl
julia
```
and activate the package by running
```
]activate .
```
You should now be able to run the tests and examples /test and /examples
respectively, and the complete documentation can be found in an interactive
HTML-format [here](http://www.control.lth.se/personnel/marcus-greiff/), as generated from docs/DOCS.md. Running
```
include("examples/example_FSF_attitude.jl")
```
produces the following output

![Lyapunov function](/docs/images/robust_SU2_FSF_analysis.png "alt text")
