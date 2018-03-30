clc
clear

r_RT_inner = 0.076200; % Inner radius in meters (3 inches)
sig_y = [271 118 260]; % Yield Strength in MPa
sig_y_names = ['6061-T6      ' '6061-T6_Welded     ' '6061-T6_PWHT     ']; % Corresponding material names
n = 3;
P_atm = 0.1015; % atm pressure in MPa
P = 5.06 - P_atm; % Pressure in MPa

for i=1:3
    t(i) = (P*r_RT_inner*n)/(sig_y(i)); % designing for yield
end

t_inches = t .* 39.3701; % converting to inches

L = 0.85; % Length of cylindrical section (in meters)
A = (pi*L.*((2*r_RT_inner) + (2/1000).*t)) + (pi.*((2*r_RT_inner)+((2/1000).*t)).^2);
%RhoPERarea = [3.934 9.942 4.417]; % density per unit area in kg/m^2
%m = RhoPERarea.*A;

disp('Materials Considered: ')
disp(sig_y_names)
disp(' ')
disp('Required R.T. wall thickness (in m)')
disp(t)
disp('Required R.T. wall thickness (in inches)')
disp(t_inches)
%disp('Run tank mass (excluding weld strengtheners) in kg')
%disp(m)

