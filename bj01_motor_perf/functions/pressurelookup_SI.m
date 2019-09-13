function [pressure] = pressurelookup_SI (alt)
P = [101300, 89880, 79500, 70120, 61660, 54050, 47220, 41110, 35650, 30800, 26500, 12110, 5529];
altitude = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000];
pressure = interp1(altitude, P, alt);
end