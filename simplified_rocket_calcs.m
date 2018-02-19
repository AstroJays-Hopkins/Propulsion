clc,close all

% Assumptions
% -Thrust is constant and steady throughout burn
% -Drag coefficient of 0.75

% Constants
air_density=1.22; % density of air(kg/m^3)
drag_coef=0.75; % drag coefficient of rocket body [axial] (unitless)
D=0.075; % x-sectional diameter of rocket body [axial] (m)
g=9.81; % acceleration due to gravity (m/s^2)

% Inputs
T=input('Motor Thrust (N)? '); % Average Thrust (Newtons)
I=input('Total Impulse (N-s)? '); % Total Impulse of Rocket(Newton-seconds)

m_struc=input('Mass of Rocket Structure (kg)? '); % Mass of rocket structure
m_elec=input('Mass of Electronics/Avionics System(s) (kg)? '); % Mass of electronics/avionics system(s)
m_recov=input('Mass of Recovery System(s) (kg)? '); % Mass of recovery system (kg)
m_prop_dry=input('Mass of Propulsion excluding propellant (kg)? '); % Mass of propulsion system (dry weight)
m_propellant=input('Mass of Propellant (kg)? '); % Mass of propellant (fuel) 

%----Calculations----%
%  Rocket Properties
A=(pi()*(D^2))/4; % x-sectional area of rocket body [axial] (m^2)
m_avg_boost=(m_struc)+(m_recov)+(m_elec)+(m_prop_dry)-(m_propellant/2); % average mass during boost (kg)
m_avg_coast=(m_struc)+(m_elec)+(m_prop_dry)-(m_propellant); % average mass during coast (kg)

% Flight Parameters
t=I/T; % calculating effective burn time (seconds)
k=(0.5)*air_density*drag_coef*A; % calculating "k" (kg/m)
q=sqrt((T-(m_avg_boost*g))/k); % calculating "q" (kg/s)
q_a=sqrt((m_avg_coast*g)/k); % calculating "q_a"
q_b=sqrt((g*k)/m_avg_coast); % calculating "q_b"
x=(2*k*q)/m_avg_boost; % calculating "x" (kg/m-s)

% Parameters of Interest
v_burnout=q*((1-exp(-x*t))/(1+exp(-x*t))); % calculating burnout velocity (m/s)
y_burnout=(-m_avg_boost/(2*k))*(log((T-(m_avg_boost*g)-(k*(v_burnout^2)))/(T-(m_avg_boost*g)))); % calculating burnout altitude (m)

y_coast=(m_avg_coast/(2*k))*log(((m_avg_coast*g)+(k*(v_burnout^2)))/(m_avg_coast*g)); % calculating coasting altitude gain
t_coast=(atan(v_burnout/q_a))/q_b; % calculating coasting time

y_apogee=y_burnout+y_coast; % calculating rocket apogee (m)
y_apogee_ft=y_apogee*(3.28084); % converting apogee altitude from meters to ft


% Displaying Parameters of Interest
disp(y_apogee)
disp(y_apogee_ft)

