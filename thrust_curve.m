%% Combustion flow parameters
%CHAMBER, all parameters defined by ProPep3
P_c = 3.7e6; %Chamber pressure, Pa
T_c = 3280.713; %Chamber temp, K (from ProPep3)
k_c = 1.246; %ratio of specific heats of combustion chamber constituents, from ProPep3
k_e = 1.257; %ratio of specific heats of exhaust constituents, from ProPep3

MM = 25.781; %Molecular mass of combustion constituents, g/mol
R_prime = 8.314; % Universal gas constant,J/mol*K
R = R_prime*1000/MM; %combustion gas constant, J/kg*K

%THROAT
T_star = T_c*(1+(k_c-1)/2)^-1;
P_star = P_c*(T_star/T_c)^(k_c/(k_c-1));
c_star = sqrt(k_e*R*T_star); %speed of sound at throat
u_star = c_star; %M=1 at throat
rho_star = k_c*P_star/(c_star^2);

%EXHAUST
P_e = pressurelookup_SI(2286); %pressure at nozzle exit (chosen parameter, P(7500ft, mid-burn))
T_e = T_star*(P_e/P_star)^((k_e-1)/k_e); %NEED TO DEFINE AMBIENT PRESSURE
c_e = sqrt(k_e*R*T_e); %speed of sound at nozzle exit
u_e = sqrt(((2*k_e*R*T_star*(1-(P_e/P_star)^((k_e-1)/k_e)))/(k_e-1))+u_star^2); %gas velocity of exhaust
M_e = u_e/c_e; %Mach number of exhaust

%% Design
ER = (1/M_e)*sqrt(((1+((k_e-1)/2)*M_e^2)/(1+(k_e-1)/2))^((k_e+1)/(k_e-1))); %Expansion ratio of nozzle
OF = 7; %7:1 oxidizer:fuel ratio

t_burn = 13;
Thrust_max = 3000;
mdot_estimate = Thrust_max/u_e; %total mass flow rate estimate
m_total = mdot_estimate*t_burn;
m_fuel = m_total*(1/(OF+1));
m_ox = m_total*(OF/(OF+1));

A_star = mdot_estimate/(rho_star*u_star);
A_e = A_star*ER;
D_star = sqrt(4*A_star/3.1415);
D_e = sqrt(4*A_e/3.1415);

rdot_estimate = 0.0012; %average regression rate = 1.2 mm/s, from Aspire Space
G_max = 600; %maximum estimated mass flux through port to avoid flameout, kg/m^2*s
rin_fuel = sqrt((1/G_max)/3.1415);
rout_fuel = rin_fuel+rdot_estimate*t_burn; %outer radius of fuel, 3in

rho_fuel = 935; %kg/m^3
L_fuel = m_fuel/(rho_fuel*3.1415*((rout_fuel^2)-(rin_fuel^2)));

r_tank = 0.0762; %3in to m
U = 0.15; %percent ullage
rho_n2o_l = 772.25; %density of n2o liquid, kg/m^3
rho_n2o_v = 1.834; %density of n2o vapour, kg/m^3
L_tank = (m_ox - (rho_n2o_v*U + rho_n2o_l*(1-U))*((4/3)*3.1415*r_tank^3))/((rho_n2o_v*U + rho_n2o_l*(1-U))*(3.1415*r_tank^2));

%% thrust curve estimation, sea level
deltat = 0.1; %time step for simulation

%vector initialization
T = zeros(1, t_burn/deltat);
mdot_total = zeros(1, t_burn/deltat);
mdot_ox = m_ox/t_burn;
mdot_fuel = zeros(1, t_burn/deltat);
OFR = zeros(1, t_burn/deltat);
r = zeros(1, t_burn/deltat);
r(1) = rin_fuel;
rdot = zeros(1, t_burn/deltat);
Go = zeros(1, t_burn/deltat);

%regression rate formula: rdot = a*Go^n
a = 0.002265;
n = 1;

correction_factor = 0.975; %correction factor for thrust knockdown on 15degree half angle conical nozzle vs ideal bell nozzle

for t = 2:t_burn/deltat
   Go(t) = mdot_ox/(3.1415*r(t-1)^2); %calculate oxizider mass flux
   rdot(t) = (a*Go(t)^n)/1000;
   
   mdot_fuel(t) = 970*2*3.1415*r(t-1)*rdot(t)*L_fuel;
   mdot_total(t) = mdot_ox + mdot_fuel(t);
   T(t) = (mdot_total(t) * u_e)*correction_factor - (P_e - 101000) * A_e; %Calculate thrust
   r(t) = r(t-1) + rdot(t)*deltat; %Calculate new radius of fuel port
   
   OFR(t) = mdot_ox/mdot_fuel(t); %Oxidizer to fuel ratio
end

I_total = 0;
for i = 1:t_burn/deltat
    I_total = I_total + T(i)*deltat;
end
    
Isp = T(2)/(mdot_total(2)*9.81)

figure(1)
hold on
plot((1:t_burn/deltat)*deltat, T, 'r');
axis([0, 9, 0, max(T)]);
xlabel('Time (s)');
ylabel('Thrust (N)')
hold off

figure(2)
plot((1:t_burn/deltat)*deltat, OFR, 'r');
axis([0, 9, 0, max(OFR)]);
xlabel('Time (s)');
ylabel('Oxidizer to Fuel Ratio')
hold off

figure(3)
plot((1:t_burn/deltat)*deltat, rdot*1000, 'r');
axis([0, 9, 0, max(rdot)*1000]);
xlabel('Time (s)');
ylabel('Regression rate (mm/s)')
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
m_structure = 15; %mass of structure and engine
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
            F_T(t) = T(t) - (P_e - 101000) * A_e + (P_e - pressurelookup_SI(altitude_FS(t-1))) * A_e;
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

figure(4)
hold on
plot(time, altitude_FS, 'r')
axis([0, flight_time, 0, max(altitude_FS)]);
xlabel('time (s)');
ylabel('altitude (m)');
hold off

g_force = acceleration_FS/9.81;

figure(5)
hold on
plot(time, g_force, 'k');
axis([0, t_burn, min(g_force), max(g_force)]);
xlabel('time (s)');
ylabel('G Force (g)');
hold off