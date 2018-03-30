d_i = 0.14732; % Inner Diameter in meters (5.8 inches)
sig_y = [271 118 260]; % Yield Strength in MPa
sig_y_names = ['6061-T6      ' '6061-T6_Welded     ' '6061-T6_PWHT     ']; % Corresponding material names
n = 1.5;
P = 5.06; % Pressure in MPa

for i=1:3
    t(i) = sqrt(3)* (P*(d_i/2))/(2*n.*sig_y(i));
end
t_inches = t ./ 25.4; % converting to inches

%a = 3*(d_i^2);
%b = 6*(d_i);
%c = (4 - (16*(sig_y.^2)/((n^2)*(P^2))));


%u = ((b^2) + sqrt((b^2) - (4*a*c)))/(2*a);

%t = 1000./u; % required thickness in mm
%t_inches = t ./ 25.4; % converting to inches

%L = 0.85; % Length of cylindrical section (in meters)
%A = (pi*L.*(d_i + (2/1000).*t)) + (pi.*(d_i+((2/1000).*t)).^2);
%Rho_area = [3.934 9.942 4.417]; % density per unit area in kg/m^2
%m = Rho_area.*A;

disp(' ')
disp(sig_y_names)
disp(' ')
disp('Required R.T. wall thickness (in mm)')
disp(t)
disp('Required R.T. wall thickness (in inches)')
disp(t_inches)
%disp('Run tank mass (excluding weld strengtheners) in kg')
%disp(m)

