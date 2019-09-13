function [density] = densitylookup_n2o_SI(p)  %density of saturated nitrous oxide liquid

%%!!THIS IS LEGACY SOFTWARE (as of 09/13/2019) USEFUL ONLY FOR ORIGINAL MOTOR
%%SIZING SCRIPT. HENCEFORTH USE refpropm.m or N2Olookup.m for
%%thermophysical properties of N2O!!!
rho = [907.4, 881.6, 853.5, 822.2, 786.6, 743.9, 688, 589.4, 452]; % [kg/m3]
press = [3127, 3547, 4007, 4510, 5060, 5660, 6315, 7033, 7251]*1000; % [Pa]
density = interp1(press, rho, p);
end