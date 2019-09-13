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

% ----------------
% NOMENCLATURE:
% - "RT" = "Run Tank", the flight oxidizer tank
% - "liq" = refers to the bulk liquid phase of N2O in the RT
% - "ull" = refers to the "ullage" in the RT
% - "CC" = "Combustion Chamber", the 


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
Conv.MtoFT = 3.28084; %[m to ft]
Conv.MtoIN = 39.3701; %[m to in]
Conv.KGtoLBM = 2.20462; %[kg to lbm]

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
InjLine.length1 = /Conv.MtoFT; % length of first segment of vent line (the vertical portion) [ft --> m]
InjLine.length2 = /Conv.MtoFT; % length of 2nd segment of vent line (the horizontal portion) [ft --> m]
InjLine.length3 = /Conv.MtoFT; % length of outlet to amb segment of vent line (the horizontal portion) [ft --> m]

% ------ Injector Geometry & Flow Coefficients ------ %
Injector.CdA = /(Conv.MtoFt^2); % effective discharge area of injector [ft^2 --> m^2]

% ------ Combustion Chamber Geometry ------ %
Fuel.OuterR = /Conv.MtoIN; % outer radius of the fuel [in --> m]

% ------ Nozzle Geometry ------ %
Nozzle.DivAngle = deg2rad(15); % half angle of divergence of conical nozzle [rad]

%% 2. Initial Conditions (& those that are engineer variable)

% ------ Ambient Init Conditions ------ %
Amb.alt = /Conv.MtoFT; % altitude above sea-level at start of burn [ft --> m]

% ------ Run Tank Init Conditions ------ %
RT.bulk.mass = /Conv.KGtoLBM; % initial total mass of N2O in RT (both liquid and ullage) [lbm --> kg]
RT.liq.temperature = ; % initial temperature of bulk N2O liquid [°F]

% ------ Fuel Init Conditions ------ %
Fuel.InnerR = /conv.MtoIN; % initial inner radius of solid fuel port [in --> m]

% ------ Nozzle Init Conditions ------ %
Nozzle.throat.dia = ; % initial throat diameter of nozzle [ft]
Nozzle.exit.dia = ; % initial exit diameter of nozzle [ft]




