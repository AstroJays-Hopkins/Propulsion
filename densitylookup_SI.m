function [density] = densitylookup_SI(alt)
rho = [1.125, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.5900, 0.5258, 0.4671, 0.4135, 0.1948];
altitude = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000];
density = interp1(altitude, rho, alt);
end