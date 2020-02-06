function [val_desired] = N2O_NonSat_Lookup(temp,press,desired_prop)
% Takes N2O temp and pressure as well as the desired property output var as inputs
% and returns that desired output at the given TD-state

% INPUT TEMP = Kelvin
% INPUT PRESS = MPa

% desired_prop:
    % OPTIONS
        % Density = kg/m^3
        % IntEnergy = kJ/kg
        % Enthalpy = kJ/kg
        % Entropy = kJ/kg-K
        % Cv = kJ/kg-K
        % Cp = kJ/kg-K
        % CpCv = unitless (Cp/Cv)
        % SoundSpeed = m/s
        % CompFactor = unitless
        % Phase = unitless
        % ThermCond = mW/m-K
        % Viscosity = microPascals-sec
        % ThermDiff = cm^2/sec
        % Prandtl = unitless

        
Temps = [250,260,270,280,285,290,295,300,305,310];
temp_diff = Temps - temp;
[~,tempidx]=sort(abs(temp_diff));
lb_tempind = tempidx(1);
lb_temp = Temps(lb_tempind);
ub_tempind= tempidx(2);
ub_temp = Temps(ub_tempind);

Press = 0.1:0.05:12;
press_diff = Press - press;
[~,pressidx] = sort(abs(press_diff));
lb_pressind = pressidx(2);
ub_pressind= pressidx(1);

lb = strcat('N2O_',num2str(lb_temp),'K');
ub = strcat('N2O_',num2str(ub_temp),'K');

low = load('N2O_NonSat_Properties',lb);
upp = load('N2O_NonSat_Properties',ub);
low = struct2table(low);
upp = struct2table(upp);

load('N2O_NonSat_Properties','varnames');
varind = find(varnames == desired_prop);

tempinterp = [lb_temp,ub_temp];
dv(1) = interp1(Press,table2array(low.(lb)(:,varind)),press);
dv(2) = interp1(Press,table2array(upp.(ub)(:,varind)),press);

val_desired = interp1(tempinterp,dv,temp);

end