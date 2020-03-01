%% Documentation
% Homogenous Equilibrium Model (HEM) based Injector Calcs
% Forcing injector to be cavitationally choked acc. to HEM

%% Setup
clc, clear, close all
dQ = 0.0001;
s_tol = 0.0001;

% Inputs
Cd = 0.65;
A = (0.1:0.05:1)*(0.0254^2); % inj area [m^2]
T1 = FtoC(70); % inlet temp [C]

% Getting rest of upstream conditions
P1 = N2Olookup('temperature',T1,1,'pressure'); % getting pressure [kPa] (assuming saturated, pure liq) 
h1 = N2Olookup('temperature',T1,1,'enthalpy'); % getting enthalpy [kJ/kg] (assuming saturated, pure liq) 
s1 = N2Olookup('temperature',T1,1,'entropy'); % getting entropy [kJ/kg] (assuming saturated, pure liq) 

% Downstream 
P2 = linspace(0.5*P1,P1,30); % chamber pressure [kPa]
s2 = s1; % isentropic process acc. to methodology of HEM

for p = 1:length(P2)
    for i = 1:length(A)
        conv = false;
        if i == 1
            Q_guess = 0.3;
        else
            Q_guess = Q(p,i-1);
        end
        ii(p,i) = 1;
        while conv == false
            s_guess = N2Olookup('pressure',P2(p),Q_guess,'entropy');
            err = s1 - s_guess;
            if err < -s_tol
                Q_guess = Q_guess - dQ;
                ii(i) = ii(p,i) + 1;
            elseif err > s_tol
                Q_guess = Q_guess + dQ;
                ii(i) = ii(p,i) + 1;
            else
                Q(p,i) = Q_guess;
                conv = true;
            end
        end
        rho2(p,i) = N2Olookup('pressure',P2(p),Q(p,i),'density');
        h2(p,i) = N2Olookup('pressure',P2(p),Q(i),'enthalpy');
        mdot(p,i) = Cd*A(i)*rho2(i)*sqrt(2*1000*(h1-h2(i))); %[kg/s]
    end
end
%%
pr = P2/P1;
plot(pr,mdot(:,5))