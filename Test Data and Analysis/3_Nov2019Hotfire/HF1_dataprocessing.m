%% Documentation

% Data Processing from Hotfire 1 ("HF1") of BJ-01 Hybrid Engine
% Original Author: Dan Zanko (11/22/2019)

%% Initialization

clear, clc, close all

ft = load('HF1_fulltest.mat');
fb = load('HF1_fullburn.mat');
lp = load('HF1_nomburn.mat');


% NOTE::::: 
%      structure tag "ft" = full test data
%      structure tag "fb" = full burn data (includes long trailoff)
%      structure tag "lp" = liquid phase portion of burn before blowdown

% ** this was done for ease of coding/shorthand **

%% Engine Parameters

% Fuel Grain
    eng.rho_fuel = 0.033779; % density of HDPE [lbm/in^3]
    eng.L = 38.64; % total initial fuel grain length [in.]
    eng.r_in = 1; % initial inner radius of fuel [in.]
    eng.r_out = 1.5; % outer radius of fuel [in.]

% Injector
    eng.d_injhole = 5/64; % [in]
    eng.n_injholes = 14;
    eng.A_injhole = pi*(eng.d_injhole^2)/4; % [in^2]
    eng.A_inj = eng.n_injholes * eng.A_injhole; % [in^2]


%% Creating important vars from data to use in calcs later on

% Key Points in Burn 
    fb.liqdepl = 106; % index corresponding to depletion of liquid in RT
    fb.end = 158; % index corresponding to thrust hitting ~10% of the eyeballed avg. val during liq. phase of burn

% var for easy referencing of length
    ft.l = length(ft.LCR1);
    ft.dtime = diff(ft.time);

    fb.l = length(fb.LCR1);
    fb.dtime(1,1) = 0;
    fb.dtime(2:fb.l,1) = diff(fb.time);

    lp.l = length(lp.LCR1);
    lp.dtime(1,1) = 0;
    lp.dtime(2:lp.l,1) = diff(lp.time);
    
% N2O Mass Flow Rates
    fb.mdot_ox(1,1) = 0;
    fb.mdot_ox(2:fb.l,1) = -diff(fb.LCR1)./fb.dtime(2:fb.l,1);

    lp.mdot_ox(1,1) = 0;
    lp.mdot_ox(2:lp.l,1) = -diff(lp.LCR1)./lp.dtime(2:lp.l,1);
    
% Thrust (total thrust = sum of the two thrust load cells)
    fb.thrust = fb.LCC1 + fb.LCC2;
    lp.thrust = lp.LCC1 + lp.LCC2;

%% Ox Flow Calcs

% WARNING THAT THIS IS EFFECTIVELY TAKING THE CdA of the ENTIRE LOWER
% PLUMBING LINES (NOT the injector)

fb.deltaP = fb.PTR1 - fb.PTC1; % diff in pressure from ox tank to CC (psi)

for i = 1:fb.l
    fb.rho_ox(i,1) = 0.000036127*N2O_NonSat_Lookup(FtoK(fb.TCR3(i)),psi2MPa(fb.PTR1(i)),"Density"); % density of n2o downstream of runtank converted to lbm/in^3
    if fb.deltaP(i,1) < 0
        fb.CdA(i,1) = 0;
    else
        fb.CdA(i,1) = fb.mdot_ox(i,1)./sqrt(2*32.2*12*fb.rho_ox(i,1).*fb.deltaP(i,1));
    end
end


%% Engine Calcs

% Chamber Pressure
fb.PchAvg.liq = mean(fb.PTC1(1:fb.liqdepl));
fb.PchAvg.ull = mean(fb.PTC1(fb.liqdepl:fb.end));
fb.PchAvg.tot = mean(fb.PTC1(1:fb.end));

% Oxidizer Spent
fb.SpentOx.liq = fb.LCR1(1) - fb.LCR1(fb.liqdepl);
fb.SpentOx.ull = fb.LCR1(fb.liqdepl) - fb.LCR1(fb.end);
fb.SpentOx.tot = fb.SpentOx.liq + fb.SpentOx.ull;

% Ox. Mass Flow Rate ((m_fin - m_init) / deltat)
fb.AvgMdotOx.liq = (fb.LCR1(1) - fb.LCR1(fb.liqdepl)) / fb.time(fb.liqdepl);
fb.AvgMdotOx.ull = (fb.LCR1(fb.liqdepl) - fb.LCR1(fb.end)) / (fb.time(fb.end) - fb.time(fb.liqdepl));
fb.AvgMdotOx.tot = (fb.LCR1(1) - fb.LCR1(fb.end)) /f b.time(fb.end);

% Thrust
fb.maxthrust.liq = max(fb.thrust(1:fb.liqdepl));
fb.maxthrust.ull = max(fb.thrust(fb.liqdepl:fb.end));
fb.maxthrust.tot = max(fb.thrust(1:fb.end));

fb.thrustavg.liq = mean(fb.thrust(1:fb.liqdepl));
fb.thrustavg.ull = mean(fb.thrust(fb.liqdepl:fb.end));
fb.thrustavg.tot = mean(fb.thrust(1:fb.end));

% Impulse
fb.totimpulse.liq = sum(fb.dtime(1:fb.liqdepl).*fb.thrust(1:fb.liqdepl));
fb.totimpulse.ull = sum(fb.dtime(fb.liqdepl:fb.end).*fb.thrust(fb.liqdepl:fb.end));
fb.totimpulse.tot = fb.totimpulse.liq+fb.totimpulse.ull;


%% (SEE DISCLAIMER) Regression Calcs

% **DISCLAIMER** ---------------
% Because we actually ran out of HDPE to combust, these regression calcs
% are incredibly sketchy

% Ammount of fuel spent per phase (liq and vap injections) was calculated
% by taking the total fuel spent and then saying that the fuel spent in
% each of these two phases is weighted by the fraction of the total impulse
% developed at this phase of the burn

fb.fuelspent.tot = 4.93; % fuel spent over total burn (inital - measured final val) [lbm]
fb.impulsefrac = fb.totimpulse.liq/fb.totimpulse.tot; % fraction of the total impulse that was generated during the liq phase [-]
fb.fuelspent.liq = fb.impulsefrac*fb.fuelspent.tot;
fb.fuelspent.ull = fb.fuelspent.tot - fb.fuelspent.liq;

fb.r_ld = sqrt( (fb.fuelspent.liq / (eng.rho_fuel*pi*eng.L)) + (eng.r_in^2) ); % getting fuel inner rad @ liq depletion [in.]
fb.rdot_avg.liq = (eng.r_in - fb.r_ld) / fb.time(fb.liqdepl); % average rdot over course of entire liq phase of burn [in/s]

% From mass_ox vs. time plot we see we have an approx. linear rate of
% decrease, allowing us to simplify to approx. constant mdot_ox during the
% liq phase of the burn



% **END OF DISCLAIMER RELATED CALCS**

%% -------------------------------- Plots ------------------------------------- %%
    %% Thrust vs. time
    figure('Name', 'Thrust'), hold on

        plot(fb.time,fb.thrust, fb.time,smooth(fb.thrust));

        xlabel('Time (s)')
        ylabel('Thrust (lbf)')
        title('Thrust vs. Time - HF1')
        grid on, grid minor
        legend('Raw','Smoothed')

    hold off  
    %% Ox Mass and tank pressure vs. time
    figure, hold on
        subplot(2,1,1)
            plot(fb.time,fb.LCR1)
            grid on, grid minor
            xlabel('Time (s)')
            ylabel('N_2O Mass (lbm)')
            title('N_2O Mass vs. Time - HF1')
        subplot(2,1,2)
            plot(fb.time,fb.PTR1)
            grid on, grid minor
            xlabel('Time (s)')
            ylabel('Tank Pressure (psia)')
            title('Lower Tank Pressure vs. Time - HF1')
    hold off
    %% Plotting thrust and chamber pressure vs. time
    figure('Name', 'Thrust & Chamber P'), hold on
        
        xlabel('Time (s)')
        title('Thrust and Chamber Pressure vs. Time - HF1')
        grid on, grid minor
        
        yyaxis left
            plot(fb.time,fb.thrust);
            ylabel('Thrust (lbf)')
        yyaxis right
            plot(fb.time,fb.PTC1); 
            ylabel('Chamber Pressure (psia)')
            ylim([0 inf])

    hold off
    %% Plotting ox tank and chamber pressures
    figure('Name', 'Tank & Chamber P'), hold on

        plot(fb.time, fb.PTR1, fb.time, fb.PTC1);

        xlabel('Time (s)')
        ylabel('Pressure (psia)')
        title('Pressure vs. Time - HF1')
        legend('Tank Pressure','Chamber Pressure')
        grid on, grid minor

    hold off   
    %% Plotting mdot and deltaP b/w OxTank and CC vs. time 
    figure('Name', 'mdot & dP'), hold on
        subplot(2,1,1)
            plot(fb.time,fb.mdot_ox,fb.time,smooth(fb.mdot_ox))
            xlabel('Time (s)')
            ylabel('Oxidizer Mass Flow Rate (lbm/s)')
            grid on, grid minor
            title('Oxidizer Flow Rate vs. Time - HF1')
            legend('Raw','Smoothed')
            ylim([0 inf])
        subplot(2,1,2)
            plot(fb.time, fb.deltaP,fb.time,smooth(fb.deltaP));
            xlabel('Time (s)')
            ylabel('\DeltaP (psig)')
            grid on, grid minor
            title('Lower Tank-to-Chamber \Delta P vs. Time - HF1')
            legend('Raw','Smoothed')
    hold off
    %% Plotting mdot and Pch & thrust vs. time 
    figure('Name', 'mdot & dP'), hold on
        subplot(3,1,1)
            plot(fb.time,fb.mdot_ox)
            xlabel('Time (s)')
            ylabel('Oxidizer Mass Flow Rate (lbm/s)')
            grid on, grid minor
            title('(Smoothed) Oxidizer Flow Rate vs. Time - HF1')
            ylim([0 inf])
        subplot(3,1,2)
            plot(fb.time, fb.PTC1);
            xlabel('Time (s)')
            ylabel('Chamber Pressure (psia)')
            grid on, grid minor
            title('(Smoothed) Chamber Pressure vs. Time - HF1')
            legend('Raw')
        subplot(3,1,3)
            plot(fb.time, fb.thrust);
            xlabel('Time (s)')
            ylabel('Thrust (lbf)')
            grid on, grid minor
            title('(Smoothed) Thrust vs. Time - HF1')
            legend('Raw')
    hold off 
    %% Ullage and Liq Temps and tank pressure vs. time
    figure, hold on
        subplot(2,1,1)
            plot(fb.time,fb.TCR3, fb.time,fb.TCR1)
            grid on, grid minor
            xlabel('Time (s)')
            ylabel('Temperature (ºF)')
            title('N_2O Run Tank Temperatures vs. Time - HF1')
            legend('Tank Outlet Fluid Temp', 'Ullage Temp')
        subplot(2,1,2)
            plot(fb.time,fb.PTR1)
            grid on, grid minor
            xlabel('Time (s)')
            ylabel('Tank Pressure (psia)')
            title('Lower Tank Pressure vs. Time - HF1')
    hold off
    %% Combustion Chamber Temps and Chamber pressure vs. time
    figure, hold on
        subplot(2,1,1)
            plot(fb.time,fb.TCC1,fb.time,fb.TCC2,fb.time,fb.TCC3)
            grid on, grid minor
            xlabel('Time (s)')
            ylabel('Temperature (ºF)')
            title('Combustion Chamber Temperatures vs. Time - HF1')
            legend('Injector','Pre-Comb Chamber','Post-Comb Chamber')
        subplot(2,1,2)
            plot(fb.time,fb.PTC1)
            grid on, grid minor
            xlabel('Time (s)')
            ylabel('Tank Pressure (psia)')
            title('Lower Tank Pressure vs. Time - HF1')
    hold off
    %% (CURRENTLY COMMENTED OUT) plotting effective discharge area vs. time
    % 
    % 
    % figure('Name', 'CdA'), hold on
    % 
    %     plot(fb.time,lp.CdA)
    %     
    %     fb.Cd = fb.CdA/eng.A_inj;
    %     xlabel('Time (s)')
    %     ylabel('Effective Discharge Area (in^2)')
    %     title('Effective Discharge Area vs. Time - HF1 - TAKEN USING P_{TANK} b/c INJECT PT FAILED - USE WITH CAUTION')
    %     grid on, grid minor
    %     
    % hold off