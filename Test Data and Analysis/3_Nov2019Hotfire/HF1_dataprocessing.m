%% Documentation % Initialization
% Data Processing from Hotfire 1 ("HF1") of BJ-01 Hybrid Engine
% Original Author: Dan Zanko (11/22/2019)

clear, clc, close all

ft = load('HF1_fulltest.mat');
nb = load('HF1_nomburn.mat');

% NOTE::::: 
%      structure "ft" = full test
%      structure "nb" = nominal burn
% ** this was done for ease of coding/shorthand **

%% Creating important vars from data to use in calcs later on

% Ask max for data entry for "ignite signal send"
% Ask max for "MV-G1 moving" data entry
ft.l = length(ft.LCR1); % var for easy referencing of length
ft.thrust = ft.LCC1 + ft.LCC2; % total thrust is the sum of the load cell readings
ft.dtime = diff(ft.time);

nb.l = length(nb.LCR1); % var for easy referencing of length
nb.dtime(1,1) = 0;
nb.dtime(2:nb.l,1) = diff(nb.time);
nb.thrust = nb.LCC1 + nb.LCC2; % total thrust is the sum of the load cell readings
nb.mdot_ox(1,1) = 0;
nb.mdot_ox(2:nb.l,1) = -diff(nb.LCR1)./nb.dtime(2:nb.l,1);
nb.mdot_ox(7,1) = 0; % taking out a weird highly neg datapoint

%% Engine Performance Calcs

nb.totimpulse = sum(nb.dtime.*nb.thrust);
nb.totimpulse_extrap = nb.totimpulse*36/30; %extrapolating what the total impulse would be for a burn using a full tank
nb.thrustavg = mean(nb.thrust(7:nb.l));

%% Ox Flow Calcs
d_injhole = 5/64; % [in]
n_injholes = 14;
A_injhole = pi*(d_injhole^2)/4; % [in^2]
A_inj = n_injholes*A_injhole; % [in^2]

nb.deltaP = nb.PTR1 - nb.PTC1; % diff in pressure from ox tank to CC (psi)
for i = 1:nb.l
    nb.rho_ox(i,1) = 0.000036127*N2O_NonSat_Lookup(FtoK(nb.TCR3(i)),psi2MPa(nb.PTR1(i)),"Density"); % density of n2o downstream of runtank converted to lbm/in^3
    if nb.deltaP(i,1) < 0
        nb.CdA(i,1) = 0;
    else
        nb.CdA(i,1) = nb.mdot_ox(i,1)./sqrt(2*32.2*12*nb.rho_ox(i,1).*nb.deltaP(i,1));
    end
end


%% Plots

%% plotting thrust vs. time
figure('Name', 'Thrust'), hold on

    plot(nb.time, nb.thrust);
    
    xlabel('Time (s)')
    ylabel('Thrust (lbf)')
    title('Thrust vs. Time - HF1')
    grid on, grid minor
    
hold off

%% plotting thrust and chamber pressure vs. time
figure('Name', 'Thrust & Chamber P'), hold on
    
    yyaxis left
        plot(nb.time,nb.thrust);
        ylabel('Thrust (lbf)')
    
    yyaxis right
        plot(nb.time,nb.PTC1); 
        ylabel('Chamber Pressure (psia)')
    
    xlabel('Time (s)')
    title('Thrust and Chamber Pressure vs. Time - HF1')
    grid on, grid minor
    
hold off

%% plotting ox tank and chamber pressures
figure('Name', 'Tank & Chamber P'), hold on
    
    plot(nb.time, nb.PTR1, nb.time, nb.PTC1);
    
    xlabel('Time (s)')
    ylabel('Pressure (psia)')
    title('Pressure vs. Time - HF1')
    legend('Tank Pressure','Chamber Pressure')
    grid on, grid minor
    
hold off

%% plotting mdot and deltaP b/w OxTank and CC vs. time 
figure('Name', 'mdot & dP'), hold on

    yyaxis left
    	plot(nb.time, nb.mdot_ox)
        ylabel('Oxidizer Mass Flow Rate (lbm/s)')
        
    yyaxis right
        plot(nb.time, nb.deltaP);
        ylabel('Pressure Drop Across Main Ox Line & Injector (psia)')
        
    xlabel('Time (s)')
    title('Ox Flow Rate and \Delta P vs. Time - HF1')
    legend('mdot_{ox}','\Delta P_{tank-combch}')
    grid on, grid minor
    
hold off

%% plotting effective discharge area vs. time
figure('Name', 'CdA'), hold on

    plot(nb.time,nb.CdA)
    
    nb.Cd = nb.CdA/A_inj;
    xlabel('Time (s)')
    ylabel('Effective Discharge Area (in^2)')
    title('Effective Discharge Area vs. Time - HF1')
    grid on, grid minor
    
hold off
