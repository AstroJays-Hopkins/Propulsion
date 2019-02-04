function [Cd] = CDcurve (Ma)

Ma_points = [0.0, 0.3, 0.6, 0.8, 1.0, 1.3, 1.5];
Cd_points = [0.25, 0.21, 0.25, 0.3, 0.6, 0.6, 0.6];

if Ma > 1.5
    Ma = 1.5;
end
    
Cd = interp1(Ma_points, Cd_points, Ma);
end