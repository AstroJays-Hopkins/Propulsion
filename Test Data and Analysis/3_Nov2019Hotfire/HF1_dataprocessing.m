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
nb.dtime = diff(nb.time);
nb.thrust = nb.LCC1 + nb.LCC2; % total thrust is the sum of the load cell readings
nb.mdot_ox = diff(nb.LCR1)./nb.dtime;

%% 
nb.totimpulse = sum(nb.dtime.*nb.thrust(2:nb.l));
nb.totimpulse_extrap = nb.totimpulse*36/30; %extrapolating what the total impulse would be for a burn using a full tank
nb.thrustavg = mean(nb.thrust(7:nb.l));
nb.deltaP = nb.PTR1 - nb.PTC1; % diff in pressure from ox tank to CC (psi)

%% Plots

% plotting thrust vs. time
figure, hold on
    thrtime = plot(nb.time, nb.thrust);
    xlabel('Time (s)')
    ylabel('Thrust (lbf)')
    title('Thrust vs. Time - HF1')
    grid on, grid minor
hold off

% plotting thrust and chamber pressure vs. time
figure, hold on
    thrpchtime = plotyy(nb.time,nb.thrust,nb.time,nb.PTC1); 
    xlabel('Time (s)')
    ylabel(thrpchtime(1),'Thrust (lbf)')
    ylabel(thrpchtime(2),'Chamber Pressure (psia)')
    title('Thrust and Chamber Pressure vs. Time - HF1')
    grid on, grid minor
hold off

% plotting ox tank and chamber pressures
figure, hold on
    presstime = plot(nb.time, nb.PTR1, nb.time, nb.PTC1);
    xlabel('Time (s)')
    ylabel('Pressure (psia)')
    title('Pressure vs. Time - HF1')
    legend('Tank Pressure','Chamber Pressure')
    grid on, grid minor
hold off

% plotting mdot and deltaP b/w OxTank and CC vs. time 
figure, hold on
    mdotdptime = plotyy(nb.time(2:nb.l), nb.mdot_ox, nb.time(2:nb.l), nb.deltaP(2:nb.l));
    xlabel('Time (s)')
    ylabel(mdotdptime(1),'Oxidizer Mass Flow Rate (lbm/s)')
    ylabel(mdotdptime(2),'Pressure Drop Across Main Ox Line & Injector (psia)')
    title('Ox Flow Rate and deltaP vs. Time - HF1')
    legend('mdot_{ox}','\Delta P_{tank-combch}')
    grid on, grid minor
hold off
