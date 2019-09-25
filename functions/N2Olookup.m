function [value] = N2Olookup(input_type,input,quality,prop)

%% Documentation
% Uses NIST table for saturated (Liq-Vap) N2O to get a desired property
% ASSUMES SATURATED PURE N2O

% Inputs:
%       - input_type = the type of parameter used for input
            % OPTIONS
            % 'temperature' = temperature [°C]
            % 'pressure'    = pressure [kPa]
            
%       - input = value of the input prop

%       - quality = defines which state you want the property for
            % 0 = pure liquid
            % 1 = pure vapor
                 
%       - prop = the property you want
            % OPTIONS
            % 'temperature' = Saturation temperature [°C] (INPUT TYPE MUST BE PRESSURE)
            % 'pressure' .  = Saturation pressure [kPa] (INPUT TYPE MUST BE TEMPERATURE)
            % 'cp'          = Const Pressure Specific Heat Capacity [kJ/kg*K]
            % 'cv'          = Const Volume Specific Heat Capacity [kJ/kg*K]
            % 'density'     = Density [kg/m^3]
            % 'enthalpy'    = mass specific enthalpy [kJ/kg]
            % 'entropy'     = mass specific entropy [kJ/kg*K]
            % 'JouleThom'   = !!NOTE UNITS!! Joule-Thomson expansion coefficient [F/psia]
            % 'pressure'    = Saturation Pressure [kPa]
            % 'SoundSpd'    = speed of sound [m/s]
            % 'SpVol'       = mass-specific volume [m^3/kg]
            % 'SurfTen'     = Surface tension coefficient [N/m]

%% -------- MAIN CODE -------- %%.
load 'N2O_SatProperties.mat'; % Loading dataset

if strcmp(prop,'temperature') % checks if requesting saturation temperature
    if strcmp(input_type, 'pressure') % checks to make sure a pressure was input
        value = interp1(Sat.pressure, Sat.temperature, input); % gets sat temp
    else
        fprintf('FAILURE, to get saturation temperature the input type must be a saturation pressure')
    end
    
elseif  strcmp(prop,'pressure') == 1 % checks if requesting saturation pressure
    if strcmp(input_type, 'temperature') % checks to make sure a temperature was input
        value = interp1(Sat.temperature, Sat.pressure, input); % gets sat press
    else
        fprintf('FAILURE, to get saturation pressure the input type must be a saturation temperature')
    end
    
else % in the case the user is requesting neither saturation pressure nor temperature
    if strcmp(input_type,'temperature') == 1
        liq_value = interp1(Sat.temperature, Liq.(prop), input);
        vap_value = interp1(Sat.temperature, Vap.(prop), input);
    elseif strcmp(input_type, 'pressure') == 1
        liq_value = interp1(Sat.pressure, Liq.(prop), input);
        vap_value = interp1(Sat.pressure, Vap.(prop), input);
    end
    value = ((1-quality)*liq_value) + (quality*vap_value);
end

end