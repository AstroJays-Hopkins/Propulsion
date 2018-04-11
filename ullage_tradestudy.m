clc, clear, close all

temp_step = 1;
T_close = (0:temp_step:25); %RT temperature upon final closing of bleed valve [deg C]
T_close = transpose(T_close);
T_maxop = 25; %Max designed for RT temp [deg C]
U_maxop = 0.01; % setting volumetric percent of vapor in RT at T_maxop
[rho_RT_liq_maxop, rho_RT_vap_maxop, P_RT_vap_maxop] = N2Olookup(T_maxop);

LENGTH = length(T_close);
rho_RT_liq_close = zeros(LENGTH,1);
rho_RT_vap_close = zeros(LENGTH,1);
P_RT_vap_close = zeros(LENGTH,1);
U = zeros(LENGTH,1);

i=1;

while i <= LENGTH
    [rho_RT_liq_close(i), rho_RT_vap_close(i), P_RT_vap_close(i)] = N2Olookup(T_close(i));
    U(i) = ( ( ((1 - U_maxop) * rho_RT_liq_maxop ) + (U_maxop*rho_RT_vap_maxop) - rho_RT_liq_close(i) ) / (rho_RT_vap_close(i) - rho_RT_liq_close(i)) ); % percent ullage
    i = i +1;
end

figure, hold on
plot(T_close,U), xlabel('Temperature at Tank Closing (Celcius)'), ylabel('Ullage Fraction Required')

% vol_tank = m_ox./((rho_RT_liq_close.*(1-U))+(rho_RT_vap_close.*U)); % calculates required volume of tank (m^3)
% r_tank = 3*(2.54/100); % setting inner run tank radius to 3" (converting to meters)
% L_tank_hemiends = (vol_tank-((4/3)*pi*(r_tank^3)))/(pi*(r_tank^2)); % assumes cylindrical tank with 2 hemispherical end caps
% L_tank_flatends = vol_tank/(pi*(r_tank^2)); % assumes cylindrical tank with 2 flat end caps (w/ no volume contribution)