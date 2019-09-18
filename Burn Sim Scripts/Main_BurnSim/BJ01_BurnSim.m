%% Documentation
% Engine burn simulation script
% BJ-01 Hybrid Engine Development - HDPE + N2O
% Written by Dan Zanko - 09/12/2019
% Last Updated: 09/12/2019

% ----------------------
% STRUCTURE OF SCRIPT:
% 0. Unit conversion constants
% 1. Inputting Engine Geometry (at start of burn)
% 2. Initial Conditions of Burn
% 3. Simulation Options/Config
% 4. Initialization of Simulation

% ----------------
% NOMENCLATURE:
% - "RT" = "Run Tank", the flight oxidizer tank
% - "liq" = refers to the bulk liquid phase of N2O in the RT
% - "ull" = refers to the "ullage" in the RT
% - "CC" = "Combustion Chamber", the assembly which houses the solid fuel
% - "amb" = ambient conditions
% - "traj" = flight trajectory

% -------------
% NOTES:
% - Solid fuel configuration is that of a single, axially aligned
%   cylindrical port
% - Nozzle is conical


clear, clc, close all
addpath /BJ-01 Motor Sims/Burn Sim Scripts
addpath /BJ-01 Motor Sims/functions
addpath /N2O Properties

%% 0. Unit Conversions
Conv.PSItoPa = 6894.74; %[psi to Pa]
Conv.LBFtoN = 4.44822; %[lbf to N]
Conv.MtoFt = 3.28084; %[m to ft]
Conv.MtoIN = 39.3701; %[m to in]
Conv.KGtoLbm = 2.20462; %[kg to lbm]

%% 1. Inputting Engine Geometry (NON-variable system parameters)

% ------ Loss Coefficients ------ %
% All of these are "K" values for the component/config
loss.teeStraight = ; % loss coeff assoc. with flow through a tee going straight through
loss.tee90 = ; % loss coeff assoc. with flow through a bent line, changing dir by 90°
loss.bend = ; % loss coeff assoc. with flow through a bent line, changing dir by 90°
loss.mNPTtoTube = ; % loss coeff
loss.NOSolenoid = ; % loss coeff of N.O. solenoid used in sys
loss.NCSolenoid = ; % loss coeff of N.C. solenoid used in sys {NOS Super Pro Shot}

% ------ Vent Line Geometry & Flow/Loss Coefficients ------ %
VentLine.length1 = /Conv.MtoFT; % length of first segment of vent line (the vertical portion) [ft --> m]
VentLine.length2 = /Conv.MtoFT; % length of 2nd segment of vent line (the horizontal portion) [ft --> m]
VentLine.length3 = /Conv.MtoFT; % length of outlet to amb segment of vent line (the horizontal portion) [ft --> m]

% ------ Run Tank Geometry ------ %
RT.bulk.vol = /(Conv.MtoFt^3); % total internal volume of N2O flight tank [ft^3 --> m^3]
RT.bulk.height = /Conv.MtoFT; % total height of internal volume of RT [ft --> m]

% ------ Injection Line Geometry & Flow/Loss Coefficients ------ %
% InjLine.length1 = /Conv.MtoFT; % length of first segment of vent line (the vertical portion) [ft --> m]
% InjLine.length2 = /Conv.MtoFT; % length of 2nd segment of vent line (the horizontal portion) [ft --> m]
% InjLine.length3 = /Conv.MtoFT; % length of outlet to amb segment of vent line (the horizontal portion) [ft --> m]
InjLine.CdA = 0.000129511324641318/(Conv.MtoFt^2); % effective discharge area of injector line [ft^2 --> m^2]

% ------ Injector Geometry & Flow Coefficients ------ %


% ------ Combustion Chamber Geometry ------ %
Fuel.OuterR = /Conv.MtoIN; % outer radius of the fuel [in --> m]

% ------ Nozzle Geometry ------ %
Nozzle.DivAngle = deg2rad(15); % half angle of divergence of conical nozzle [rad]

%% 2. Initial Conditions (& those that are engineer variable)

% ------ Ambient Init Conditions ------ %
amb.alt = 152/Conv.MtoFT; % altitude above sea-level at start of burn [ft --> m]

% ------ Run Tank Init Conditions ------ %
RT.bulk.mass = /Conv.KGtoLBM; % initial total mass of N2O in RT (both liquid and ullage) [lbm --> kg]
RT.bulk.temperature = FtoK( ); % initial temperature of bulk N2O liquid [°F-->K]

% ------ Fuel Init Conditions ------ %
Fuel.InnerR = /Conv.MtoIN; % initial inner radius of solid fuel port [in --> m]

% ------ Nozzle Init Conditions ------ %
Nozzle.throat.dia = 0.966/Conv.MtoIN; % throat diameter [in-->m]
Nozzle.exit.dia = 2.171/Conv.MtoIN; % exit diameter [in-->m]

%% 3. Simulation Config 

tstep = 0.1; % delta-t we use for our time-stepping simulation
flight = false; % setting if a static hotfire or a flight sim ("false" and "true" respectively) 

%% 4. Initialization of Sim
fuel = true; % creating a boolean to track if there's still fuel in rocket
ox = true;  % creating boolean to track if there's still 
i = 1; % creating iterator

while(fuel == true && ox == true) % simulation runs as long as there's both fuel and oxidizer left in the engine
    
    amb.pressure(i) = pressurelookup_SI(amb.alt(i)); % getting ambient pressure at current altitude [Pa]
    % Liquid state in RT
    RT.pressure(i) = N2Olookup("temperature",RT.bulk.temperature-273.15,0,"pressure")*1000; % assuming saturated liquid, getting pressure [kPa-->Pa]
    RT.rho(i) = N2Olookup("temperature",RT.bulk.temperature-273.15,0,"density"); % assuming saturated liquid, getting density [kg/m^3]
    
    % chamber pressure guess
    if i == 1
        CC.pressguess(i) = 400*Conv.PSItoPa; % initial chamber pressure guess of 400 psi [psi-->Pa]
    else
        CC.pressguess(i) = CC.pressure(i-1);

    end

    % mass flow rate & flux
    CC.mdot.ox(i) = InjLine.CdA * sqrt( 2*RT.rho(i)*(RT.pressure(i)-CC.pressguess(i)) ); % mass flow rate of oxidizer into combustion chamber [kg/s]
    
    
end

