################################################################################
# Copyright (C) Marcus Greiff 2020
#
# @file example_attitude_sys_id.jl
# @author Marcus Greiff
# @date July 2020
#
# This scrips contains some data and methods used in the system ideantification
# of Patrik's home-made drone. The basic methods and ideas are quite genetral,
# but not documeneted here.
#
################################################################################
using LinearAlgebra
using LaTeXStrings

# In the first experiment, we measure the thrust of the rotors as a function of
# a normalized command sent to the motors (typically pwm duty cycle), here,
# pairs of (measured force, duty cycle) are given in F and D respectively, from
# which the a and b parameters of the polynomial map are using a LS-regression.

# Generate some synthetic data
Theta    = rand(2,1);
N        = 10^3;
D        = rand(N,1);
Phi      = [D D.^2];
F        =  Phi*Theta + 0.1*randn(N,1);


#length = 55.0;

D = hcat(collect(50:50:600))
W = hcat([
    6,
    19,
    37,
    60,
    87,
    119,
    148,
    180,
    211,
    242,
    271,
    305,
])

g = 9.81;
m = 560;

D = D./1000;
W = W./1000;

#f = ma

F = W.*g * 2 / 4.0;


Phi      = [D D.^2];
Cov      = inv(Phi'*Phi);
Thetahat = (Phi'*Phi)\(Phi'*F);

using Plots
dd = hcat(collect((0:N)/N))
scatter(D, F, marker=:circle, color=:green, markeralpha = 0.3, xaxis="Duty cycle", yaxis="Force", label="Measurements")
#plot!(dd, [dd dd.^2]*Theta, color=:red, linewidth=2, label="True function")
plot!(dd, [dd dd.^2]*Thetahat, color=:blue, linewidth=2, label="Estimated function")

#1.1885640614704445
#2.291161604353094


# In the first experiment, we estimate C using the above function and
# measurements of the rotational torque about the e3 axis as a function of the
# applied PWM duty cycle. Here, we whould expect a linear relationship between
# the rotor thrust and the generated torque. If this is not linear, then c
# should represent the linear approximation of this function at the single rotor
# thrust of mg/4.


tauz = 4*0.005*9.81*0.09; # Nm
fi   = 2*1.18701; # N
# tauz = c * fi where fi is the thrust of one of the rotors
c = tauz / fi;

# index 0 - front right - pos force, positive moment about z, negative moment x, negative moment y
# index 1 - back right  - pos force, negative moment about z, negative moment x, positive moment y
# index 2 - back left   - pos force, positive moment about z, positive moment x, positive moment y
# index 3 - front left  - pos force, negative moment about z, positive moment x, negative moment y
#
#     3       0
#       .   .
#         .
#       .   .
#     2       1
#
# long 15.5
# short 9.3

A = [+1.0, +1.0, +1.0, +1.0;
     -dy, -dy, +dy, +dy;
     -dx, +dx, +dx, -dx;
     + c, - c, + c, - c]

# Cylindrical cross (mass, length of one axis, radius of the cross)
MC = 0.284
LC = 0.24
RC = 0.01

ICX = MC * (RC^2 / 4 + LC^2 / 12) + MC*RC^2 / 2;
ICY = ICX;
ICZ = MC * (RC^2 / 4 + LC^2 / 12) + MC * (RC^2 / 4 + LC^2 / 12);

# Battery box (Mass, length (x), with (y), height (z), offset in z)
ME = 0.172
WE = 0.035
LE = 0.09
HE = 0.025
DE = 0.03

IEX = ME * (WE^2/12 + HE^2/12 + DE^2);
IEY = ME * (LE^2/12 + HE^2/12 + DE^2);
IEZ = ME * (LE^2/12 + HE^2/12);

IXX = ICX + IEX;
IYY = ICY + IEY;
IZZ = ICZ + IEZ;
