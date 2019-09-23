function [temp_K] = FtoK(temp_F)
% Converts an input temperature from °F to Kelvin

temp_K = ((temp_F-32)*(5/9))+273.15;

end

