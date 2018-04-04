function [rho_RT_liq, rho_RT_vap, P_RT_vap] = N2Olookup(temp)

% Defining datasets
RT_temperature = [0,5,10,15,20,25,30,35,36.42]; % run tank temperature [deg. C]
N2O_liq_rho = [907.4, 881.6, 853.5, 822.2, 786.6, 743.9, 688.0, 589.4, 452.0]; % nitrous oxide liquid density [kg/m^3]
N2O_vap_rho = [86.7, 100.2, 116.1, 135.4, 159.4, 191.1, 237.3, 330.5, 452.0]; % nitrous oxide vapor density [kg/m^3]
N2O_vap_pressure = [3.13, 3.55, 4.01, 4.51, 5.06, 5.66, 6.31, 7.03, 7.25]; % nitrous oxide vapor pressure [MPa]

% Interpolation
rho_RT_liq = interp1(RT_temperature, N2O_liq_rho, temp);
rho_RT_vap = interp1(RT_temperature, N2O_vap_rho, temp);
P_RT_vap = interp1(RT_temperature, N2O_vap_pressure, temp);

end