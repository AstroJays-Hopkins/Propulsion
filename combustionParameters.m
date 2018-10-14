% Function that accepts the current oxidizer-to-fuel ratio and Chamber
% Pressure, and outputs the temperature in the chamber [K], the ratio of
% specific heats, and the material-specific gas constant. Only use for
% combustion of Nitrous Oxide and HDPE

function [T, gamma, R] = combustionParameters (OF, P_c)

% Chamber Pressures corresponding to each row of properties, [psi]
P_c_points = [350, 525];

%
OF_points = [5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10];

% Chamber Temperature, [K]
T_points = [2894, 3040, 3144, 3211, 3249, 3267, 3273, 3270, 3261, 3249, 3235;
            2905, 3058, 3168, 3242, 3285, 3306, 3313, 3310, 3302, 3289, 3274];
   
% Ratio of Combustion Products Specific Heats
gamma_points = [1.2622, 1.2562, 1.2518, 1.2487, 1.2467, 1.2454, 1.2446, 1.2441, 1.2439, 1.2439, 1.2439;
                1.2618, 1.2556, 1.2510, 1.2477, 1.2455, 1.2441, 1.2432, 1.2428, 1.2425, 1.2425, 1.2426];

% Molar Mass, kg/mol
MM_points = [23.283, 24.050, 24.704, 25.254, 25.714, 26.101, 26.430, 26.715, 26.962, 27.179, 27.371;
             23.304, 24.084, 24.753, 25.317, 25.787, 26.182, 26.517, 26.080, 27.052, 27.269, 27.460];
         
R_points = 8314./MM_points; % [J/kg K]


if P_c > 525
    P_c = 525; % Update tables with more representative parameter space in future
elseif P_c < 350
    P_c = 350;
end

if OF > 10
    OF = 10; % Update tables with more representative parameter space in future
elseif OF < 5
    OF = 5;
end
        
T1 = interp1(OF_points, T_points(1,:), OF);
T2 = interp1(OF_points, T_points(2,:), OF);
T = interp1(P_c_points, [T1, T2], P_c);

gamma1 = interp1(OF_points, gamma_points(1,:), OF);
gamma2 = interp1(OF_points, gamma_points(2,:), OF);
gamma = interp1(P_c_points, [gamma1, gamma2], P_c);

R1 = interp1(OF_points, R_points(1,:), OF);
R2 = interp1(OF_points, R_points(2,:), OF);
R = interp1(P_c_points, [R1, R2], P_c);

end