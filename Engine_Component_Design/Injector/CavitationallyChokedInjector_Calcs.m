%% Documentation
% Homogenous Equilibrium Model (HEM) based Injector Calcs
% Forcing injector to be cavitationally choked acc. to HEM

%% Setup
% Inputs
Cd = 0.65;
A = (0.1:0.01:0.35)*(0.0254^2); % inj area [m^2]
T1 = FtoC(70); % inlet temp [C]

% Getting rest of upstream conditions
P1 = N2Olookup('temperature',T1,1,'pressure'); % getting pressure [kPa] (assuming saturated, pure liq) 
h1 = N2Olookup('temperature',T1,1,'enthalpy'); % getting enthalpy [kJ/kg] (assuming saturated, pure liq) 
s1 = N2Olookup('temperature',T1,1,'entropy'); % getting entropy [kJ/kg] (assuming saturated, pure liq) 

% Downstream
P2 = psi2MPa(500)/1000; % chamber pressure [MPa]
s2 = s1; % isentropic process acc. to methodology of HEM

%%
% for i = 1:length(A)
% %    rho2 = N2O_NonSat_Lookup(T1,P2,
%    mdot(i) = Cd*A(i)*rho2*sqrt(2*(h1-h2(i))); 
% end
% mdot_crit = max(mdot);
% plot(A, mdot)