%BJ-01 Hybrid Engine Development - HDPE + N2O
%written by Andrew Colombo
%COPYRIGHT ASTROJAYS 2018
%
%Structure of script:
%1. Initialization of constants and vectors
%2. Initial isentropic calculations to estimate engine size parameters
%3. Engine size estimates
%4. Thrust curve estimate
%5. Recalculate engine sizing based on actual fuel requirements
%6. Calculate required run tank sizing
%7. Plot engine performance
%8. Flight sim
%9. Plot flight simulation results

%% INITIALIZATION
%constant initialization
C_characteristic = 1550; %[m/s] characteristic exhaust velocity (set by propep3)
k_c = 1.246; %[] ratio of specific heats of combustion chamber constituents, from ProPep3
k_e = 1.257; %[] ratio of specific heats of exhaust constituents, from ProPep3
P_drop = 0.6666666; %[] fraction of Run Tank pressure after emptying
t_burn = 13; %[s]
deltat = 0.01; %[s] time step for simulation


MM = 25.781; %[g/mol] Molecular mass of combustion constituents
R_prime = 8.314; % [J/mol*K] Universal gas constant
R = R_prime*1000/MM; %[J/kg*K] combustion gas constant

%vector creation
T = zeros(1, t_burn/deltat); %[N] THRUST
mdot_total = zeros(1, t_burn/deltat); %[kg/s] MASS FLOW, TOTAL
mdot_fuel = zeros(1, t_burn/deltat); %[kg/s] MASS FLOW, FUEL
mdot_ox = zeros(1, t_burn/deltat); %[kg/s] MASS FLOW, OX
OFR = zeros(1, t_burn/deltat); %[] OXIDIZER TO FUEL RATIO
r = zeros(1, t_burn/deltat); %[m] PORT RADIUS 
rdot = zeros(1, t_burn/deltat); %[m/s] REGRESSION RATE
Go = zeros(1, t_burn/deltat); %[kg/s] OXIDIZER MASS FLUX
P_c = zeros(1, t_burn/deltat); %[Pa] CHAMBER PRESSURE
P_RT = zeros(1, t_burn/deltat); %[Pa] RUN TANK PRESSURE
u_e = zeros(1, t_burn/deltat); %[m/s] EXHAUST EXIT VELOCITY
P_e_TC = zeros(1, t_burn/deltat); %[Pa] EXHAUST PRESSURE
deltaP = zeros(1, t_burn/deltat); %[Pa] PRESSURE DIFFERENT BETWEEN RT AND CHAMBER



%% Combustion flow CONSTANTS & Initial Isentropic Calculations to Size Engine
%Initialize chamber pressure at target value and calculate isentropic flow
%conditions down nozzle to get estimates to size engine

%CHAMBER, all parameters defined by ProPep3
P_c(1) = 3.7e6; %[Pa] Initialize chamber pressure at target value, Pa
T_c = (3291+3240)/2; %[K] Average Chamber temp (from ProPep3) Note: chamber temp does not vary significantly with changing OF ratio and chamber pressure

%THROAT
T_star = T_c*(1+(k_c-1)/2)^-1;
P_star = P_c(1)*(T_star/T_c)^(k_c/(k_c-1));
c_star = sqrt(k_e*R*T_star); %speed of sound at throat
u_star = c_star; %M=1 at throat
rho_star = k_c*P_star/(c_star^2);

%EXHAUST
P_e = pressurelookup_SI(1500); %pressure at nozzle exit (P(1500m, mid-burn))
T_e = T_star*(P_e/P_star)^((k_e-1)/k_e); %NEED TO DEFINE AMBIENT PRESSURE
c_e = sqrt(k_e*R*T_e); %[m/s] speed of sound at nozzle exit
u_e(1) = sqrt(((2*k_e*R*T_star*(1-(P_e/P_star)^((k_e-1)/k_e)))/(k_e-1))+u_star^2); %gas velocity of exhaust
M_e = u_e(1)/c_e; %Mach number of exhaust

%% Design of motor based on estimated performance parameters calculated above
ER = (1/M_e)*sqrt(((1+((k_e-1)/2)*M_e^2)/(1+(k_e-1)/2))^((k_e+1)/(k_e-1))); %Expansion ratio of nozzle, expanding to pressure at estimated altitude for mid-burn
OF = 7; %7:1, initial oxidizer to fuel ratio

Thrust_max = 3000;
mdot_estimate = Thrust_max/u_e(1); %total mass flow rate estimate
m_total = mdot_estimate*t_burn;
m_fuel = m_total*(1/(OF+1));
m_ox = m_total*(OF/(OF+1));

A_star = mdot_estimate/(rho_star*u_star); %[m2]
A_e = A_star*ER; 
D_star = sqrt(4*A_star/3.1415); %[m]
D_e = sqrt(4*A_e/3.1415);

rdot_estimate = 0.0012; %[m/s] average regression rate = 1.2 mm/s, from Aspire Space
G_max = 600; %[kg/m^2*s] maximum estimated mass flux through port to avoid flameout
rin_fuel = sqrt((1/G_max)/3.1415);
rout_fuel = rin_fuel+rdot_estimate*t_burn; %[m] outer radius of fuel, 3in

rho_fuel = 935; %[kg/m^3] Density of HDPE, 
L_fuel = m_fuel/(rho_fuel*3.1415*((rout_fuel^2)-(rin_fuel^2)));


rho_n2o_l = 786.6; %FOR T=20ºC, density of n2o liquid, kg/m^3
rho_n2o_v = 159.4; %FOR T=20ºC, density of n2o vapour, kg/m^3

CR = 4; %Contraction ratio of nozzle

%% thrust curve estimation, sea level
%Vector initialization
r(1) = rin_fuel;
P_c(1) = 0.1013e6; %Combustion chamber pressure equal to atmospheric at start
P_RT(1) = 5.06e6; %Run tank pressure, 50.6 bar
deltaP(1) = P_RT(1) - P_c(1);

%Injector orifice diameter calc
r_inj = sqrt((m_ox/t_burn)/(sqrt(2*rho_n2o_l*((5.06-3.7+0)*10^6))*pi)); %[m] calculates injector orifice radius
A_inj = pi()*(r_inj)^2; %[m^2] Injector area

%Estimated constants for regression rate formula: rdot = a*Go^n
a = 0.002265;
n = 1;

%Knockdown factors
Conical_Nozzle_Correction_Factor = 0.983; %correction factor for thrust knockdown on 15degree half angle conical nozzle vs ideal bell nozzle
Chamber_Throat_Area_Ratio_Knockdown = 0.99; %reduction in thrust due to losses in converging section of nozzle

%Emptying of Run Tank
P_drop_per_step = (P_RT(1)-P_RT(1)*P_drop)/(t_burn/deltat); %Pressure drop of the Run Tank per time step

for t = 1:((t_burn/deltat)-1)
   mdot_ox(t) = A_inj*sqrt(2*(deltaP(t))*rho_n2o_l);
   
   Go(t) = mdot_ox(t)/(pi*r(t)^2); %calculate oxizider mass flux
   rdot(t) = (a*Go(t)^n)/1000;
   
   mdot_fuel(t) = rho_fuel*2*pi*r(t)*rdot(t)*L_fuel;
   mdot_total(t) = mdot_ox(t) + mdot_fuel(t);
   
   
   %ISENTROPIC CALCS:
   %THROAT
   T_star_TC = T_c*(1+(k_c-1)/2)^-1;
   P_star_TC = P_c(t)*(T_star_TC/T_c)^(k_c/(k_c-1));
   c_star_TC = sqrt(k_e*R*T_star_TC); %speed of sound at throat
   u_star_TC = c_star_TC; %M=1 at throat
   rho_star_TC = k_c*P_star_TC/(c_star_TC^2);
   %EXHAUST
   P_e_TC(t) = pressurelookup_SI(1500); %pressure at nozzle exit (P(1500m, mid-burn))
   T_e_TC = T_star*(P_e_TC(t)/P_star_TC)^((k_e-1)/k_e); %NEED TO DEFINE AMBIENT PRESSURE
   c_e_TC = sqrt(k_e*R*T_e_TC); %speed of sound at nozzle exit
   u_e(t) = sqrt(((2*k_e*R*T_star_TC*(1-(P_e_TC(t)/P_star_TC)^((k_e-1)/k_e)))/(k_e-1))+u_star_TC^2); %gas velocity of exhaust
   
   
   T(t) = (mdot_total(t) * u_e(t))*Conical_Nozzle_Correction_Factor*Chamber_Throat_Area_Ratio_Knockdown - (P_e_TC(t) - 101300) * A_e; %Calculate thrust
   
   r(t+1) = r(t) + rdot(t)*deltat; %Calculate new radius of fuel port
   
   OFR(t) = mdot_ox(t)/mdot_fuel(t); %Oxidizer to fuel ratio
   
   P_RT(t+1) = P_RT(t)-P_drop_per_step;
   P_c(t+1) = mdot_total(t)*C_characteristic/A_star;
   deltaP(t+1) = P_RT(t+1)-P_c(t+1);

   %Without the conditional below, the chamber pressure oscillates too
   %wildly - resulting in the chamber pressure momentarily reaching above
   %the Run Tank pressure. This causes deltaP to be negative, which in turn
   %causes mdot_ox to have an imaginary component.
   
   %Averages out the oscillations of the chamber pressure during startup to avoid code breaking
   
   epsilon = 0.001e6; %adjust to achieve smaller oscillations
   
   if abs(P_c(t+1) - P_c(t)) > epsilon
       deltaP(t+1) = (deltaP(t)+deltaP(t+1))/2;
   end
   
end

I_total = 0;
for i = 1:t_burn/deltat
    I_total = I_total + T(i)*deltat;
end

Isp = T(2)/(mdot_total(2)*9.81); %ProPep3 

%% RECALCULATE ENGINE SIZING PARAMETERS 
%Based on usage during burn
rout_fuel = r(((t_burn/deltat)-1)); %set outer radius of propellant to the final burn radius
m_fuel = rho_fuel*L_fuel*pi*((rout_fuel^2)-(rin_fuel^2));

m_ox = 0;
for t = 1:((t_burn/deltat)-1)
    m_ox = m_ox + mdot_ox(t) * deltat; %calculate the total oxidizer used
end

%% Run Tank Sizing & "Ullage"/Head Space Calculations
ox_mass_compensation_factor = 1.25;
T_close = 20; %RT temperature upon final closing of bleed valve [deg C]
T_maxop = 25; %Max designed for RT temp [deg C]
[rho_RT_liq_close, rho_RT_vap_close, P_RT_vap_close] = N2Olookup(T_close);
[rho_RT_liq_maxop, rho_RT_vap_maxop, P_RT_vap_maxop] = N2Olookup(T_maxop);
U_maxop = 0.01; % setting volumetric percent of vapor in RT at T_maxop
U = ( ((1-U_maxop)*rho_RT_liq_maxop) + (U_maxop*rho_RT_vap_maxop) - rho_RT_liq_close) / (rho_RT_vap_close - rho_RT_liq_close); % percent ullage

r_tank = 3.125*(2.54/100); % setting inner run tank radius to 3" (converting to meters)
vol_tank = m_ox*ox_mass_compensation_factor/((rho_RT_liq_close*(1-U))+(rho_RT_vap_close*U)); % calculates required volume of tank (m^3)
vol_hemicap = (2/3)*pi*(r_tank^3);
L_tank_hemiends = (vol_tank-(2*vol_hemicap))/(pi*(r_tank^2)); % assumes cylindrical tank with 2 hemispherical end caps
h_ullage = (vol_tank*(1-U) - vol_hemicap)/(pi*r_tank^2); % calculating long. distance from BOTTOM of CYLINDER SECTION of RT the ullage valve is to be located


%% PLOT - ENGINE PERFORMANCE
figure(1)
subplot(3,1,1)
hold on
plot((1:t_burn/deltat)*deltat, T, 'r');
axis([0, 9, 0, 5000]);
title('BJ-01 Thrust vs Time');
xlabel('Time (s)');
ylabel('Thrust (N)')

subplot(3,1,2)
plot((1:t_burn/deltat)*deltat, OFR, 'r');
axis([0, 9, 0, 12]);
title('OF Ratio');
xlabel('Time (s)');
ylabel('Oxidizer to Fuel Ratio')

subplot(3,1,3)
plot((1:t_burn/deltat)*deltat, rdot*1000, 'r');
axis([0, 9, 0, 1.6]);
title('Regression Rate')
xlabel('Time (s)');
ylabel('Regression rate (mm/s)')
hold off

figure(2)
plot((1:t_burn/deltat)*deltat, P_c, 'r')
hold on
plot((1:t_burn/deltat)*deltat, P_RT, 'c')
plot((1:t_burn/deltat)*deltat, deltaP, 'k')
axis([0, 9, 0, 6e6]);
title('Run Tank & Combustion Pressure');
xlabel('Time (s)')
ylabel('Pressure (Pa)');
legend('Chamber Pressure', 'Run Tank Pressure', 'Pressure Difference b/w RT and Chamber');
hold off

figure(3)
plot((1:t_burn/deltat)*deltat, mdot_ox, 'r')
hold on
plot((1:t_burn/deltat)*deltat, mdot_fuel, 'c')
plot((1:t_burn/deltat)*deltat, mdot_total, 'k')
axis([0, 9, 0, 2]);
title('Ox and Prop Mass Flow Rates');
xlabel('Time (s)')
ylabel('Mass Flow Rate (kg/s)');
legend('Oxidizer MFR', 'Propellant MFR');
hold off


%% Flight Sim
flight_time = 45; %seconds
timesteps = flight_time/deltat;
time = linspace(0, timesteps*deltat, timesteps);

acceleration_FS = zeros(1, timesteps);
altitude_FS = zeros(1, timesteps); %initial altitude
velocity_FS = zeros(1, timesteps); %initial velocity
CA = 3.1415 * 0.0762^2; %6in DIA, cross sectional area of rocket
CD = 0.75; %from https://spaceflightsystems.grc.nasa.gov/education/rocket/shaped.html

m_pl = 4; %payload mass, kg
m_structure = 25; %mass of structure and engine
m_fuelox_FS = m_fuel + m_ox;

m_total_FS = zeros(1, timesteps);
m_total_FS(1) = m_pl + m_structure + m_fuelox_FS;

F_total = zeros(1, timesteps);
F_T = zeros(1, timesteps);
F_d = zeros(1, timesteps);

for t = 2:timesteps
    if altitude_FS(t-1) >= 0
        m_total_FS(t) = m_pl + m_structure + m_fuelox_FS;
        
        if t <= length(T)
            F_T(t) = T(t) - (P_e_TC(t) - 101000) * A_e + (P_e_TC(t) - pressurelookup_SI(altitude_FS(t-1))) * A_e;
        end

        if velocity_FS(t-1) ~= 0
            F_d(t) = 0.5 * CD * densitylookup_SI(altitude_FS(t-1)) * CA * (velocity_FS(t-1)^2)*(-velocity_FS(t-1)/abs(velocity_FS(t-1)));
        end

        F_g = -m_total_FS(t) * 9.81;
        F_total(t) = F_T(t) + F_g + F_d(t);
        acceleration_FS(t) = F_total(t)/m_total_FS(t);
        velocity_FS(t) = velocity_FS(t-1) + (acceleration_FS(t) * deltat);
        altitude_FS(t) = altitude_FS(t-1) + (velocity_FS(t) * deltat);

        if t <= length(T)
            m_fuelox_FS = m_fuelox_FS - mdot_total(t)*deltat;
        end
        
    elseif altitude_FS(t-1) < 0
        break
    end
        
end

Ma = zeros(1, timesteps);
for t = 1:timesteps
    c = sqrt(1.4 * pressurelookup_SI(altitude_FS(t))/densitylookup_SI(altitude_FS(t)));
    Ma(t) = velocity_FS(t)/c;
end

g_force = acceleration_FS/9.81;

%% PLOT - FLIGHT
figure(4)
hold on
subplot(3,1,1)
plot(time, altitude_FS, 'r')
axis([0, flight_time, 0, max(altitude_FS)]);
title('BJ-01 Flight Trajectory');
xlabel('time (s)');
ylabel('altitude (m)');

subplot(3,1,2)   
plot(time, Ma, 'r');
axis([0, flight_time, 0, 2]);
title('Mach Number');
xlabel('time (s)');
ylabel('Ma');

subplot(3,1,3)
plot(time, g_force, 'k');
axis([0, flight_time, min(g_force), max(g_force)]);
title('G Forces')
xlabel('time (s)');
ylabel('G Force (g)');
hold off