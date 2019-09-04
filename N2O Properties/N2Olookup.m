function [value] = N2Olookup(temp,quality,prop)

%% Documentation
% Uses NIST table for saturated (Liq-Vap) N2O to get a desired property
% ASSUMES SATURATED PURE N2O

% Inputs:
%       - temp = temperature of the N2O in CELSIUS

%       - quality = defines which state you want the property for
            % 0 = pure liquid
            % 1 = pure vapor
                 
%       - prop = the property you want
            % OPTIONS
            % 'cp'        = Const Pressure Specific Heat Capacity [kJ/kg*K]
            % 'cv'        = Const Volume Specific Heat Capacity [kJ/kg*K]
            % 'density'   = Density [kg/m^3]
            % 'enthalpy'  = mass specific enthalpy [kJ/kg]
            % 'entropy'   = mass specific entropy [kJ/kg*K]
            % 'JouleThom' = !!NOTE UNITS!! Joule-Thomson expansion coefficient [F/psia]
            % 'pressure'  = Saturation Pressure [kPa]
            % 'SoundSpd'  = speed of sound [m/s]
            % 'SpVol'     = mass-specific volume [m^3/kg]
            % 'SurfTen'   = Surface tension coefficient [N/m]

%% -------- MAIN CODE -------- %%.
load 'N2O_PhysicalProperties.mat'; % Loading dataset

if strcmp(prop,'pressure') == 0 % checks if requesting saturation pressure
    liq_prop = Liq.(prop);
    liq_value = interp1(Sat.Temperature, liq_prop, temp);
    vap_prop = Vap.(prop);
    vap_value = interp1(Sat.Temperature, vap_prop, temp);
    value = ((1-quality)*liq_value) + (quality*vap_value);
else % in the case the user is requesting saturation pressure
    value = interp1(Sat.Temperature, Sat.Pressure, temp);
end

end