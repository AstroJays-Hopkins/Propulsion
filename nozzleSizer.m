%NOZZLE SIZER
%written by Andrew Colombo
%COPYRIGHT ASTROJAYS 2018

load('BJ01_motor_params_061918.mat');

MtoIN = 39.3701; %[m to in]

alpha_inlet = 45; %[deg] convergence half angle of nozzle inlet
alpha_outlet = 15; %[deg] convergence half angle of nozzle outlet

L_throat = 0.0625; %[in] length of the throat section of the nozzle

L_inlet = MtoIN * (rout_fuel-D_star) / tand(alpha_inlet)

L_outlet = MtoIN * (D_e-D_star) / (2*tand(alpha_outlet))