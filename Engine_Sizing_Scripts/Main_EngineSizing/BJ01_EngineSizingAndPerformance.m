%BJ-01 Hybrid Engine Development - HDPE + N2O
%Written by Andrew Colombo & Dan Zanko
%Updated Mar 16, 2019
%
%BOTH SI & IMPERIAL UNITS
%
%Structure of script:
%0. Unit conversion constants
%1. Initialization of constants and vectors
%2. Initial isentropic calculations to estimate engine size parameters
%3. Engine size estimates
%4. Thrust curve estimate
%5. Recalculate engine sizing based on actual fuel requirements
%6. Calculate required run tank sizing
%7. Plot engine performance
%8. Flight sim
%9. Plot flight simulation results

%% 0. Unit Conversions
PSItoPa = 6894.74; %[psi to Pa]
LBFtoN = 4.44822; %[lbf to N]
MtoIN = 39.3701; %[m to in]
KGtoLBM = 2.20462; %[kg to lbm]

%% 1. Initialization 
%constant initialization, SI units due to ProPep3 output
k_c = 1.246; %[] ratio of specific heats of combustion chamber constituents, from ProPep3
k_e = 1.257; %[] ratio of specific heats of exhaust constituents, from ProPep3
P_drop = 0.6666666; %[] fraction of Run Tank pressure after emptying
t_burn = 13; %[s]
deltat = 0.01; %[s] time step for simulation

P_RT_i = 801 * PSItoPa; %starting RT pressure [Pa]
P_RT_f = P_RT_i * P_drop; %end RT pressure [Pa]


MM = 25.781; %[g/mol] Molecular mass of combustion constituents
R_prime = 8.314; % [J/mol*K] Universal gas constant
R = R_prime*1000/MM; %[J/kg*K] combustion gas constant

%vector creation
T = zeros(1, t_burn/deltat); %[N] THRUST
mdot_total = zeros(1, t_burn/deltat); %[kg/s] MASS FLOW, TOTAL
mdot_nozzle = zeros(1, t_burn/deltat); %[kg/s] MASS FLOW THROUGH NOZZLE, TOTAL
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
deltaP = zeros(1, t_burn/deltat); %[Pa] PRESSURE DIFFERENCE BETWEEN RT AND CHAMBER - Minor Losses
deltaP_minor = zeros(1, t_burn/deltat); %[Pa] Minor losses
C_characteristic = zeros(1, t_burn/deltat); %[m/s]\

%% 2. Combustion flow constants & Initial Isentropic Calculations to Grab Initial Engine Parameters
%Initialize chamber pressure at target value and calculate isentropic flow
%conditions down nozzle to get estimates to size engine

%CHAMBER CONDITIONS, all parameters defined by ProPep3
P_c(1) = 525 * PSItoPa; %[psi -> Pa] Initialize chamber pressure at target value, [Pa]
T_c = (3291+3240)/2; %[K] Average Chamber temp (from ProPep3) Note: chamber temp does not vary significantly with changing OF ratio and chamber pressure

%THROAT CONDITIONS
T_star = T_c*(1+(k_c-1)/2)^-1;
P_star = P_c(1)*(T_star/T_c)^(k_c/(k_c-1));
c_star = sqrt(k_e*R*T_star); %speed of sound at throat
u_star = c_star; %M=1 at throat
rho_star = k_c*P_star/(c_star^2);

%EXHAUST CONDITIONS
P_e = pressurelookup_SI(0); % [Pa] pressure at nozzle exit (assuming elevation of 0m)
T_e = T_star*(P_e/P_star)^((k_e-1)/k_e); %[K] temperature at nozzle exit
c_e = sqrt(k_e*R*T_e); %[m/s] speed of sound at nozzle exit
u_e(1) = sqrt(((2*k_e*R*T_star*(1-(P_e/P_star)^((k_e-1)/k_e)))/(k_e-1))+u_star^2); %[m/s] gas velocity at nozzle exit
M_e = u_e(1)/c_e; %[-] Mach number at nozzle exit

%% 3. Design of motor based on estimated performance parameters calculated above
ER = (1/M_e)*sqrt(((1+((k_e-1)/2)*M_e^2)/(1+(k_e-1)/2))^((k_e+1)/(k_e-1))); %Expansion ratio of nozzle, expanding to pressure at estimated altitude for mid-burn
OF = 7; %7:1, initial oxidizer to fuel ratio

Thrust_max = 750 * LBFtoN; %[lbf to N] set the max thrust

mdot_estimate = Thrust_max/u_e(1);
A_star_guess = mdot_estimate/(rho_star*u_star); %[m2]
A_e = A_star_guess*ER; 
D_star = sqrt(4*A_star_guess/3.1415); %[m]
D_e = sqrt(4*A_e/3.1415);

flux_SF = 1.75; %safety factor for ox mass flux
G_max = 870/flux_SF; %[kg/m^2*s] maximum estimated mass flux through port to avoid flameout with SF, based on AspireSpace literature which says 870 kg/m2.s
rin_fuel = sqrt((1/G_max)/pi);

rho_fuel = 935; %[kg/m^3] Density of HDPE

rdot_estimate = 0.0012; %[m/s] ESTIMATE average regression rate = 1.2 mm/s, from Aspire Space
rout_fuel = 1.5*0.0254; %[in -> m]setting max fuel diameter to 3in in meters
m_total = mdot_estimate*t_burn; %[kg] ESTIMATE total fuel + ox
m_fuel = m_total*(1/(OF+1)); %[kg] ESTIMATE mass of fuel
m_ox = m_total*(OF/(OF+1)); %[kg] ESTIMATE mass of oxidizer
L_fuel = m_fuel/(rho_fuel*pi*((rout_fuel^2)-(rin_fuel^2))); % [m] calculating required fuel length

% %% 4. Thrust curve estimation (at sea level)
% %Vector initialization
% A_star_prev = 0; %initialize tracker for throat area
% iter = 1; %initialize throat area iterative tracker
% A_star(1) = A_star_guess;
% 
% %Injector orifice diameter calc
% rho_n2o_l = densitylookup_n2o_SI(P_RT_i);
% C_D = 0.9;  %dischange coeff
% r_inj = sqrt((m_ox/t_burn)/(sqrt(2*rho_n2o_l*((5.06-3.7+0)*10^6))*pi)); %[m] calculates injector orifice radius
% A_inj = (pi()*(r_inj)^2)/0.9; %[m^2] Injector area
% 
% % Ballistic constants taken from Waxman dissertation
% a = 5e-2;
% n = 0.5;
% 
% %Knockdown factors
% Conical_Nozzle_Correction_Factor = 0.983; %correction factor for thrust knockdown on 15° half angle of divergence conical nozzle vs ideal bell nozzle
% Chamber_Throat_Area_Ratio_Knockdown = 0.99; %reduction in thrust due to losses in converging section of nozzle
% Combustion_Efficiency = 1;
% 
% %Emptying of Run Tank
% 
% while true
%     r = zeros(1, t_burn/deltat); %[m] PORT RADIUS 
%     r(1) = rin_fuel;
%     
%     P_c = zeros(1, t_burn/deltat); %[Pa] CHAMBER PRESSURE
%     P_c(1) = 14.6959 * PSItoPa; %[psi to Pa] Combustion chamber pressure equal to atmospheric at start
%     
%     P_RT = zeros(1, t_burn/deltat); %[Pa] RUN TANK PRESSURE
%     P_RT(1) = P_RT_i; %[Pa] Max expected run tank pressure
%     P_drop_per_step = (P_RT(1)-P_RT(1)*P_drop)/(t_burn/deltat); %Pressure drop of the Run Tank per time step
%     
%     deltaP = zeros(1, t_burn/deltat); %[Pa] Injector pressure drop
%     deltaP(1) = P_RT(1) - P_c(1);
%     
%     rho_star_TC = 2.1; %initialize for the first iteration...changes on t=2
%     C_characteristic(1) = 1550; %[m/s] characteristic exhaust velocity (set by propep3)
%     m_ox = 0;
%     
%     
%     for t = 1:((t_burn/deltat)-1)
%        rho_n2o_l = densitylookup_n2o_SI(P_RT(t));
%         
%        mdot_ox(t) = A_inj*C_D*sqrt(2*(deltaP(t))*rho_n2o_l);
%        m_ox = m_ox + mdot_ox(t) * deltat; %calculate the total oxidizer used
% 
%        Go(t) = mdot_ox(t)/(pi*r(t)^2); %calculate oxizider mass flux
%        rdot(t) = (a*Go(t)^n)/1000;
% 
%        mdot_fuel(t) = rho_fuel*2*pi*r(t)*rdot(t)*L_fuel;
%        mdot_total(t) = mdot_ox(t) + mdot_fuel(t);
% 
% 
%        %ISENTROPIC CALCS:
%        %THROAT
%        T_star_TC = T_c*(1+(k_c-1)/2)^-1;
%        P_star_TC = P_c(t)*(T_star_TC/T_c)^(k_c/(k_c-1));
%        c_star_TC = sqrt(k_e*R*T_star_TC); %speed of sound at throat
%        u_star_TC = c_star_TC; %M=1 at throat
%        if P_c < 1.0132e5
%            rho_star_TC = k_c*P_star_TC/(c_star_TC^2);
%        else
%            rho_star_TC = 2.1;
%        end
% 
%        %EXHAUST
%        P_e_TC(t) = pressurelookup_SI(0); %pressure at nozzle exit (P(1500m, mid-burn))
%        T_e_TC = T_star*(P_e_TC(t)/P_star_TC)^((k_e-1)/k_e); %NEED TO DEFINE AMBIENT PRESSURE
%        c_e_TC = sqrt(k_e*R*T_e_TC); %speed of sound at nozzle exit
%        u_e(t) = sqrt(((2*k_e*R*T_star_TC*(1-(P_e_TC(t)/P_star_TC)^((k_e-1)/k_e)))/(k_e-1))+u_star_TC^2); %gas velocity of exhaust
% 
%        T(t) = (mdot_total(t) * u_e(t))*Conical_Nozzle_Correction_Factor*Chamber_Throat_Area_Ratio_Knockdown - (P_e_TC(t) - 101300) * A_e; %Calculate thrust
% 
%        r(t+1) = r(t) + rdot(t)*deltat; %Calculate new radius of fuel port
% 
%        OFR(t) = mdot_ox(t)/mdot_fuel(t); %Oxidizer to fuel ratio
% 
%        P_RT(t+1) = P_RT(t)-P_drop_per_step;
% 
% 
%        [T_c, k_c, R] = combustionParameters(OFR(t), P_c(t));
%        if t~=1
%             C_characteristic(t) = Combustion_Efficiency*sqrt(k_c*R*T_c)/(k_c*(2/(k_c+1))^((k_c+1)/(2*k_c-2)));
%        end
% 
%        P_c(t+1) = mdot_total(t)*C_characteristic(t)/A_star(iter);
% 
%        %Minor losses
%        A_min_plumbing = 3.1415*(0.010795/2)^2;
%        v_fluid = mdot_ox(t)/(rho_n2o_l * A_min_plumbing);
%        deltaP_minor_3tree = 0.54 * rho_n2o_l * v_fluid^2 /2;
%        deltaP_minor_ball = 0.08 * rho_n2o_l * v_fluid^2 /2;
%        deltaP_minor(t) = deltaP_minor_3tree + deltaP_minor_ball;
%        
%        deltaP(t+1) = P_RT(t+1)-P_c(t+1) - deltaP_minor(t);
% 
%        %Without the conditional below, the chamber pressure oscillates too
%        %wildly - resulting in the chamber pressure momentarily reaching above
%        %the Run Tank pressure. This causes deltaP to be negative, which in turn
%        %causes mdot_ox to have an imaginary component.
% 
%        %Averages out the oscillations of the chamber pressure during startup to avoid code breaking
% 
%        epsilon = 0.001e6; %adjust to dampen pressure oscillations
% 
%        if abs(P_c(t+1) - P_c(t)) > epsilon
%            deltaP(t+1) = (deltaP(t)+deltaP(t+1))/2;
%        end
% 
%     end
%     
%     if abs(A_star_prev - A_star(iter)) < 0.0001
%         break
%     end
%         
%     A_star_prev = A_star(iter);
%     P_c_min = min(P_c(P_c>1e6));
%     OFR_min = min(OFR(OFR>0));
%     [T_c_min, k_c_min, R_min] = combustionParameters(OFR_min, P_c_min);
%     rho_c_min =  P_c_min/(R_min*T_c_min);
%     mdot_total_min = min(mdot_total(mdot_total>0));
%     u_entrance_min = mdot_total_min/(rho_c_min*pi()*rout_fuel^2);
%     Ma_entrance_min = u_entrance_min/sqrt(k_c_min*R_min*T_c_min);
%     
%     iter = iter + 1;
%     A_star(iter) = (3.1415*(rout_fuel^2)*Ma_entrance_min)*sqrt(((1+(k_c_min-1)/2)/(1+(k_c_min-1)*(Ma_entrance_min^2)/2))^((k_c_min+1)/(k_c_min-1)));
%     A_err(iter) = abs(A_star_prev-A_star(iter));
% end
% % grab min pressure & OFR in the chamber -> use to calculate the min k,
% % T_c, R -> calculate the min Ma at the entrance to the nozzle -> use
% % parameters to find the minimum throat area -> run this loop again, tracking throat area every iteration, and accept at certain convergence 
% 
% I_total = 0;
% for i = 1:t_burn/deltat
%     I_total = I_total + T(i)*deltat;
% end
% 
% Isp = T(2)/(mdot_total(2)*9.81); %ProPep3 

%% 5. Engine Sizing Parameters 
%Based on usage during burn
rout_fuel = r(((t_burn/deltat)-1)); %set outer radius of propellant to the final burn radius
m_fuel = rho_fuel*L_fuel*pi*((rout_fuel^2)-(rin_fuel^2));

%% 6. Run Tank Sizing & "Ullage"/Head Space Calculations

% SEE "BJ01_Performance_Analysis_SI.m"

%% 7. PLOT - Engine Performance
figure(1)
subplot(3,1,1)
hold on
plot((1:t_burn/deltat)*deltat, T/LBFtoN, 'r');
axis([0, 13, 0, 1000]);
title('BJ-01 Thrust vs Time');
xlabel('Time (s)');
ylabel('Thrust (lbf)')

subplot(3,1,2)
plot((1:t_burn/deltat)*deltat, OFR, 'r');
axis([0, 13, 0, 12]);
title('OF Ratio');
xlabel('Time (s)');
ylabel('Oxidizer to Fuel Ratio')

subplot(3,1,3)
plot((1:t_burn/deltat)*deltat, rdot*MtoIN, 'r');
axis([0, 13, 0, 0.07]);
title('Regression Rate')
xlabel('Time (s)');
ylabel('Regression rate (in/s)')
hold off

figure(2)
plot((1:t_burn/deltat)*deltat, P_c/PSItoPa, 'r')
hold on
plot((1:t_burn/deltat)*deltat, P_RT/PSItoPa, 'c')
plot((1:t_burn/deltat)*deltat, deltaP/PSItoPa, 'k')
plot((1:t_burn/deltat)*deltat, deltaP_minor/PSItoPa, 'k')
axis([0, 13, 0, 1000]);
title('System Pressures');
xlabel('Time (s)')
ylabel('Pressure (psi)');
legend('Chamber Pressure', 'Run Tank Pressure', 'Inj Pressure Drop', 'Minor Losses');
hold off

figure(3)
plot((1:t_burn/deltat)*deltat, mdot_ox*KGtoLBM, 'r')
hold on
plot((1:t_burn/deltat)*deltat, mdot_fuel*KGtoLBM, 'c')
plot((1:t_burn/deltat)*deltat, mdot_total*KGtoLBM, 'k')
axis([0, 13, 0, 4]);
title('Ox and Prop Mass Flow Rates');
xlabel('Time (s)')
ylabel('Mass Flow Rate (lbm/s)');
legend('Oxidizer MFR', 'Propellant MFR', 'Total MFR');
hold off

figure(4)
plot((1:t_burn/deltat)*deltat, T./(mdot_total*9.81))
hold on
axis([0, 13, 0, 300]);
title('Isp');
xlabel('Time (s)')
ylabel('Isp (s)');
hold off

%% -------------------------- %%

%% 8. Flight Sim
flight_time = 60; %seconds
timesteps = flight_time/deltat;
time = linspace(0, flight_time, timesteps);

acceleration_FS = zeros(1, timesteps);
altitude_FS = zeros(1, timesteps); %initial altitude
velocity_FS = zeros(1, timesteps); %initial velocity
CA = 3.1415 * 0.0889^2; %7in DIA, cross sectional area of rocket

m_pl = 4; %payload mass, kg
m_structure = 10.4; %mass of structure
m_propsys = 20; %dry mass of prop system
m_fuelox_FS = m_fuel + m_ox;

m_total_FS = zeros(1, timesteps);
m_total_FS(1) = m_pl + m_structure + m_propsys + m_fuelox_FS;

F_total = zeros(1, timesteps);
F_T = zeros(1, timesteps);
F_d = zeros(1, timesteps);

for t = 2:timesteps
    if altitude_FS(t-1) >= 0
        m_total_FS(t) = m_pl + m_structure + m_propsys + m_fuelox_FS;
        
        if t <= length(T)
            F_T(t) = T(t) - (P_e_TC(t) - 101000) * A_e + (P_e_TC(t) - pressurelookup_SI(altitude_FS(t-1))) * A_e;
        end
        
        
        if velocity_FS(t-1) ~= 0
            c = sqrt(1.4 * pressurelookup_SI(altitude_FS(t-1))/densitylookup_alt_SI(altitude_FS(t-1)));
            Ma = abs(velocity_FS(t-1))/c;
            
            if velocity_FS(t-1) > 0
                F_d(t) = -0.5 * CDcurve(Ma) * 1 * densitylookup_alt_SI(altitude_FS(t-1)) * CA * (abs(velocity_FS(t-1))^2);
            else
                F_d(t) = 0.5 * CDcurve(Ma) * densitylookup_alt_SI(altitude_FS(t-1)) * CA * (abs(velocity_FS(t-1))^2);
            end
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
    c = sqrt(1.4 * pressurelookup_SI(altitude_FS(t))/densitylookup_alt_SI(altitude_FS(t)));
    Ma(t) = abs(velocity_FS(t))/c;
end

g_force = acceleration_FS/9.81;

max(altitude_FS * 3.28084)

%% 9. PLOT - Flight Sim Results
figure(5)
hold on
subplot(3,1,1)
plot(time, altitude_FS*(MtoIN/12), 'r')
axis([0, flight_time, 0, 35000]);
title('BJ-01 Flight Trajectory');
xlabel('Time (s)');
ylabel('altitude (ft)');

subplot(3,1,2)   
plot(time, Ma, 'r');
axis([0, flight_time, 0, 2]);
title('Mach Number');
xlabel('Time (s)');
ylabel('Ma');

subplot(3,1,3)
plot(time, g_force, 'k');
axis([0, flight_time, min(g_force), max(g_force)]);
title('G Forces')
xlabel('Time (s)');
ylabel('G Force (g)');
hold off