%based on Aerotech N3300R-PS
clc

%% Driving Parameters
m_rocket = 13.2216; % rocket mass excluding propulsion system [kg]
m_motor_0 = 12.054; % initial motor total mass [kg]
m_prop = 7.512; % initial propellant mass [kg]

T_avg = 3168; % AVERAGE motor thrust [N]
I_tot = 14041; % total motor impulse [N-s]
t_burn = I_tot/T_avg; % burn time [s]
mdot = m_prop/t_burn;

d = 0.1016; % rocket diameter [m]
A = (pi/4)*(d^2); % rocket planform area [m^2]
C_d = 0.75; % drag coefficient of rocket [dimensionless]

rho_air = 1.22; % density of std air [kg/m^3]
g = 9.81; % accel. due to grav. [m/s^2]

%% Flight Sim
deltat = 0.01;
flight_time = 45; %seconds
timesteps = flight_time/deltat;
time = linspace(0, timesteps*deltat, timesteps);

a = zeros(1, timesteps); % acceleration
v = zeros(1, timesteps); % velocity
y = zeros(1, timesteps); % altitude

m_tot = zeros(1, timesteps);
m_tot(1) = m_rocket + m_motor_0;

F_tot = zeros(1, timesteps);
F_drag = zeros(1, timesteps);
F_W = zeros(1, timesteps);
F_t = zeros(1,timesteps);

a(1) = T_avg;
v(1) = 0;
y(1) = 0;

for t = 2:timesteps
    if t <= t_burn/deltat
        m_tot(t) = m_tot(t-1) - mdot*deltat;
        T(t) = T_avg;
    else
        m_tot(t) = m_tot(t-1);
        T(t) = 0;
    end
    
  F_drag(t) = (1/2)*C_d*A*rho_air*((v(t-1))^2);
  F_W(t) = m_tot(t)*g;
  
  F_tot(t) = T(t) - F_drag(t) - F_W(t);
  a(t) = F_tot(t)/m_tot(t);
  
  v(t) = v(t-1) + a(t-1)*deltat;
  y(t) = y(t-1) + v(t-1)*deltat;
  
end

max_accel = max(a);
max_vel = max(a);
max_alt = max(y);

y_ft = 3.28084*y; 
max_alt_ft = max(y_ft);

idx_burnout = floor(t_burn/deltat);
t_burnout = time(idx_burnout);
y_burnout = y_ft(idx_burnout);
str_burnout = ['Burnout = ', num2str(y_burnout),'ft'];

idx_apo = find(y_ft == max(y_ft));
t_apo = time(idx_apo);
y_apo = y_ft(idx_apo);
str_apo = ['Apogee = ', num2str(y_apo),'ft'];


figure
hold on
plot(time, y_ft, '-s', 'MarkerIndices', [idx_burnout idx_apo], 'MarkerFaceColor', 'red', 'MarkerSize', 15)
text(t_burnout,y_burnout,str_burnout,'HorizontalAlignment','left');
text(t_apo,y_apo,str_apo,'HorizontalAlignment','left');
xlabel('Time (s)')
ylabel('Altitude (ft)')
hold off