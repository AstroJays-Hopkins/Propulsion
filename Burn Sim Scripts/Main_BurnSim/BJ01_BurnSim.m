%% Documentation
% Engine burn simulation script
% BJ-01 Hybrid Engine Development - HDPE + N2O
% Written by Dan Zanko - 09/12/2019
% Last Updated: 09/12/2019

% -----------------------
% STRUCTURE OF SCRIPT:
% 0. Unit conversion constants
% 1. Inputting Engine Geometry (at start of burn)
% 2. Initial Conditions of Burn
% 3. Simulation Options/Config
% 4. Initialization of Simulation
% 5 Main Sim Loop

% -----------------------
% NOMENCLATURE:
% - "RT" = "Run Tank", the flight oxidizer tank
% - "liq" = refers to the bulk liquid phase of N2O in the RT
% - "ull" = refers to the "ullage" in the RT
% - "CC" = "Combustion Chamber", the assembly which houses the solid fuel
% - "amb" = ambient conditions
% - "traj" = flight trajectory
% -----------------------
% -----------------------
% OPERATION GUIDE (How to Properly Use this Script to get Engine Performance):
%
% - sim.flight = set this to:
%                       - TRUE if want to simulate a simple trajectory,
%                       - FALSE this will simply maintain initial altitude for
%                         entire burn
% -----------------------
% -----------------------
% ASSUMPTIONS:
% - Assumes that nozzle inlet conditions equivalent to stagnation
% conditions of combustion chamber
% -----------------------
% -----------------------
% NOTES:
% - Solid fuel configuration is that of a single, axially aligned
%   cylindrical port
% - Nozzle is conical
% -----------------------

%% Seting Up Workspace
clear, clc, close all
% addpath F:/Propulsion/BJ-01 Motor Sims/Burn Sim Scripts
% addpath F:/Propulsion/BJ-01 Motor Sims/functions
% addpath F:/Propulsion/N2O Properties

%% 0. Unit Conversions
Conv.PSItoPa = 6894.74; %[psi to Pa]
Conv.LBFtoN = 4.44822; %[lbf to N]
Conv.MtoFt = 3.28084; %[m to ft]
Conv.MtoIn = 39.3701; %[m to in]
Conv.KGtoLbm = 2.20462; %[kg to lbm]

%% 1. Inputting Engine Geometry (NON-variable system parameters)

% ------ Loss Coefficients ------ %
% All of these are "K" values for the component/config
loss.teeStraight = 0; % loss coeff assoc. with flow through a tee going straight through
loss.tee90 = 0; % loss coeff assoc. with flow through a bent line, changing dir by 90°
loss.bend = 0; % loss coeff assoc. with flow through a bent line, changing dir by 90°
loss.mNPTtoTube = 0; % loss coeff
loss.NOSolenoid = 0; % loss coeff of N.O. solenoid used in sys
loss.NCSolenoid = 0; % loss coeff of N.C. solenoid used in sys {NOS Super Pro Shot}

% ------ Vent Line Geometry & Flow/Loss Coefficients ------ %
VentLine.DiaOutlet = 0 / Conv.MtoIn;
VentLine.L1 = 0 / Conv.MtoFt; % length of first segment of vent line (the vertical portion) [ft --> m]
VentLine.L2 = 0 / Conv.MtoFt; % length of 2nd segment of vent line (the horizontal portion) [ft --> m]
VentLine.L3 = 0 / Conv.MtoFt; % length of outlet to amb segment of vent line (the horizontal portion) [ft --> m]

% ------ Run Tank Geometry ------ %
RT.bulk.vol = 0.751483931 /(Conv.MtoFt^3); % total internal volume of N2O flight tank [ft^3 --> m^3]
RT.bulk.height = 3.441666667 /Conv.MtoFt; % total height of internal volume of RT [ft --> m]

% ------ Injection Line Geometry & Flow/Loss Coefficients ------ %
InjLine.length1 = 0 / Conv.MtoFt; % length of first segment of vent line (the vertical portion) [ft --> m]
InjLine.length2 = 0 / Conv.MtoFt; % length of 2nd segment of vent line (the horizontal portion) [ft --> m]
InjLine.length3 = 0 / Conv.MtoFt; % length of outlet to amb segment of vent line (the horizontal portion) [ft --> m]
InjLine.CdA = 0.000129511324641318/(Conv.MtoFt^2); % effective discharge area of injector line [ft^2 --> m^2]

% ------ Injector Geometry & Flow Coefficients ------ %


% ------ Fuel Properties & Geometry ------ %
CC.fuel.OuterR = 1.5 / Conv.MtoIn; % outer radius of the fuel [in --> m]
CC.fuel.rho = 935; %[kg/m^3] Density of HDPE

% ------ Nozzle Geometry ------ %
Nozzle.DivAngle = deg2rad(15); % half angle of divergence of conical nozzle [rad]
Nozzle.ConicalConvergingKnockdown = 0.99; %reduction in thrust due to losses in converging section of nozzle
Nozzle.ConicalDivergingKnockdown = 0.983; %correction factor for thrust knockdown on 15° half angle of divergence conical nozzle vs ideal bell nozzle


%% 2. Initial Conditions & Simulation Correlating Parameters

% ------ Ambient Init Conditions ------ %
amb.alt = 152/Conv.MtoFt; % altitude above sea-level at start of burn [ft --> m]

% ------ Run Tank Init Conditions ------ %
RT.bulk.mass = 30/Conv.KGtoLbm; % initial total mass of N2O in RT (both liquid and ullage) [lbm --> kg]
RT.liq.temperature = FtoK(70); % initial temperature of bulk N2O liquid [°F-->K]

% ------ Fuel & Comb Chamber Init Conditions & Efficiencies ------ %
CC.fuel.r = 0.99 / Conv.MtoIn; % initial inner radius of solid fuel port [in --> m]
CC.fuel.L = 38.64000000 / Conv.MtoIn; % length of fuel [in --> m]
CC.cstarEfficiency = 1;
% ballistics coeffs in rdot = a(G_ox)^n model for hybrid fuel regression rate
CC.fuel.a = 5e-2; 
CC.fuel.n = 0.5;

% ------ Nozzle Init Conditions & Efficiencies ------ %
Nozzle.throat.dia = 0.966/Conv.MtoIn; % throat diameter [in-->m]
Nozzle.exit.dia = 2.171/Conv.MtoIn; % exit diameter [in-->m]

%% 3. Simulation Config 

sim.dt = 0.1; % delta-t we use for our time-stepping simulation
sim.flight = false; % setting if a static hotfire or a flight sim ("false" and "true" respectively) 
sim.P0tolerance = 100; % setting allowable deviation of chamber pressure guess from the chamber pressure calculated via nozzle theory
sim.P0change = 0.001; % if P0guess and P0 calculated deviate from each other greater than the tolerance, this value gets new guess by multiplying old guess by 1 +/- this values

%% 4. Initialization of Sim
fuel = true; % creating a boolean to track if there's still fuel in rocket
ox = true;  % creating boolean to track if there's still oxidizer in the rocket
i = 1; % creating iterator
t = 0; % creating time variable (s)

% calculating mass of loaded fuel
CC.fuel.mass = (pi/4)*(CC.fuel.OuterR^2 - CC.fuel.r^2)*CC.fuel.L*CC.fuel.rho;

% calculating areas given input diameters
Nozzle.throat.area = pi*(Nozzle.throat.dia^2)/4; % calculating nozzle throat area [m^2]

% initial oxidizer tank thermo
RT.liq.rho = N2Olookup("temperature", RT.liq.temperature-273.15, 0, "density"); % looking up density of N2O liq phase [kg/m^3]
RT.vap.rho = N2Olookup("temperature", RT.liq.temperature-273.15, 1, "density"); % looking up density of N2O vapor phase [kg/m^3]
RT.liq.volfrac = ((RT.bulk.mass/RT.bulk.vol) - RT.vap.rho)/(RT.liq.rho - RT.vap.rho); % calculating initial volume fraction of liquid phase in RT
RT.vap.volfrac = 1 - RT.liq.volfrac; % calculating initial volume fraction of vapor phase in RT
RT.liq.mass = RT.liq.volfrac*RT.bulk.vol*RT.liq.rho; % mass of liquid phase of oxidizer [kg]
RT.vap.mass = RT.liq.volfrac*RT.bulk.vol*RT.liq.rho; % mass of vapor phase of oxidizer [kg]

%% 5. Main Sim Loop
while(fuel == true && ox == true) % simulation runs as long as there's both fuel and oxidizer left in the engine
    
    amb.pressure(i) = pressurelookup_SI(amb.alt(i)); % getting ambient pressure at current altitude [Pa]
    
    % Liquid state in RT
    RT.pressure(i) = N2Olookup("temperature",RT.liq.temperature-273.15,0,"pressure")*1000; % assuming saturated liquid, getting pressure [kPa-->Pa]
    RT.liq.rho(i) = N2Olookup("temperature",RT.liq.temperature-273.15,0,"density"); % assuming saturated liquid, getting density [kg/m^3]
    
    % chamber pressure guess
    if i == 1
        CC.P0guess(i) = 400*Conv.PSItoPa; % initial chamber pressure guess of 400 psi [psi-->Pa]
    else
        CC.P0guess(i) = CC.P0(i-1);
    end
    
    % converging upon a chamber pressure
    ii = 0;
    sim.P0converge = false;
    while sim.P0converge == 0
        % mass flow rate & flux
        CC.mdot.ox(i) = InjLine.CdA * sqrt( 2*RT.liq.rho(i)*(RT.pressure(i)-CC.P0guess(i)) ); % mass flow rate of oxidizer into combustion chamber [kg/s]
        CC.oxflux(i) = CC.mdot.ox(i)/(2*CC.fuel.r(i)*pi*CC.fuel.L); % oxidizer mass flux in fuel port[kg/s-m^2]

        % fuel regression & fuel flow rate + resulting Ox/Fuel Ratio (OFR)
        CC.fuel.rdot(i) = CC.fuel.a*(CC.oxflux(i)^CC.fuel.n); % calculating regression rate of fuel [m/s]
        CC.mdot.fuel(i) = CC.fuel.rdot(i)*(2*CC.fuel.r(i)*pi*CC.fuel.L)*CC.fuel.rho; % calculating mass flow rate of fuel [kg/s]
        CC.mdot.total(i) = CC.mdot.ox(i) + CC.mdot.fuel(i);
        CC.OFR(i) = CC.mdot.ox(i)/CC.mdot.fuel(i); % calculating oxidizer to fuel ratio

        % interpolating ProPep3 outputs with current chamber pressure guess and
        % OFR to get properties of combustion products (stag temp, gamma,
        % material specific ideal gas const)
        [CC.T0(i), CC.gamma(i), CC.R(i)] = combustionParameters(CC.OFR(i), CC.P0guess(i));
        CC.char_vel(i) = (CC.cstarEfficiency/CC.gamma(i))*sqrt( (CC.gamma(i)*CC.R(i)*CC.T0(i)) / ( (2/(CC.gamma(i)+1)) ^ ((CC.gamma(i)+1)/(CC.gamma(i)-1)) ) ); % calculating characteristic velocity of exhaust
        CC.P0(i) = CC.mdot.total(i) * CC.char_vel(i) / Nozzle.throat.area(length(Nozzle.throat.area)); % calculating the chamber pressure according to nozzle theory [Pa]
        
        % checking if calculated chamber pressure matches our current guess
        % for chamber pressure
        sim.P0discrepency(i) = CC.P0(i) - CC.P0guess(i);
        if t(i) < 0.5 % enforcing P0 to remain a reasonable value on startup (prevents P0 to converge to a value >> P_tank) 
            CC.P0(i) = CC.P0guess(i);
            sim.P0converge = true;
        elseif (sim.P0discrepency(i)) < -sim.P0tolerance
            CC.P0guess(i) = CC.P0guess(i)*(1 + sim.P0change);
            ii = ii+1;
        elseif (sim.P0discrepency(i)) > sim.P0tolerance
            CC.P0guess(i) = CC.P0guess(i)*(1 - sim.P0change);
            ii = ii+1;
        else
            sim.P0converge = true; % signalling that we have converged on chamber pressure and that we can procede with the sim for this timestep
        end

    end 


    % calculating conditions at throat (choked)
    mdot_choked_book(i) = Nozzle.throat.area*CC.P0(i)*CC.gamma(i)*sqrt( (2/(CC.gamma(i)+1))^((CC.gamma(i)+1)/(CC.gamma(i)-1))) / sqrt(CC.gamma(i)*CC.R(i)*CC.T0(i));
    % calculating conditions at exhaust

    % calculating resulting thrust 

    % setting up next iteration (time-marching)
    t(i+1) = t(i) + sim.dt;
    CC.fuel.r(i+1) = CC.fuel.r(i) + CC.fuel.rdot(i)*sim.dt;
    RT.liq.mass(i+1) = RT.liq.mass(i) - CC.mdot.ox(i)*sim.dt;
    CC.fuel.mass(i+1) = CC.fuel.mass(i) - CC.mdot.fuel(i)*sim.dt;
    
    % checking if either of the propellant masses in the next time step is
    % negative, signaling to end the simulation of the burn if this
    % condition is met
    if RT.liq.mass(i+1) <= 0
        sim.ox = false;
    end
    if CC.fuel.mass(i+1) <= 0
        sim.fuel = false;
    end
    
    i = i + 1;
end

