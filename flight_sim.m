t_burn = 10;
mdot = 12.5723/t_burn; %total mass flow rate, kg/s
P_e = pressurelookup_SI(2268);
u_e = 2386; %m/s

%average thrust of MASA rocket =  2400N

timesteps = 3000; %total time = timesteps * deltat
deltat = 0.1; %time step, s
time = linspace(0, timesteps*deltat, timesteps);

a = zeros(1, timesteps);
altitude = zeros(1, timesteps); %initial altitude
velocity = zeros(1, timesteps); %initial velocity
CA = 3.1415 * 0.0762^2; %cross sectional area of rocket
CD = 0.75; %from https://spaceflightsystems.grc.nasa.gov/education/rocket/shaped.html

m_pl = 4; %payload mass, kg
m_structure = 15; %mass of structure and engine
m_fuel = mdot*t_burn; %mass of fuel for 9s burn
m_total = zeros(1, timesteps);
m_total(1) = m_pl + m_structure + m_fuel;

F_total = zeros(1, timesteps);
F_T = zeros(1, timesteps);
F_d = zeros(1, timesteps);

for t = 2:timesteps
    if altitude(t-1) >= 0
        m_total(t) = m_pl + m_structure + m_fuel;

        if m_fuel > 0
            F_T(t) = mdot*u_e+(P_e-pressurelookup_SI(altitude(t)))*6.2331*0.00053616;
        end
        
        if velocity(t-1) ~= 0
            F_d(t) = 0.5 * CD * densitylookup_SI(altitude(t-1)) * CA * (velocity(t-1)^2)*(-velocity(t-1)/abs(velocity(t-1)));
        end
        
        F_g = -m_total(t) * 9.81;
        F_total(t) = F_T(t) + F_g + F_d(t);
        a(t) = F_total(t)/m_total(t);
        velocity(t) = velocity(t-1) + (a(t) * deltat);
        altitude(t) = altitude(t-1) + (velocity(t) * deltat);

        if m_fuel > 0
            m_fuel = m_fuel - mdot*deltat;
        end
        
    elseif altitude(t-1) < 0
        break
    end
        
end

figure(1)
hold on
plot(time, altitude, 'r')
axis([0, t/10, 0, max(altitude)]);
xlabel('time (s)');
ylabel('altitude (m)');
hold off

g_force = a/9.81;

figure(2)
hold on
plot(time, g_force, 'k');
axis([0, t/10, min(g_force), max(g_force)]);
xlabel('time (s)');
ylabel('G Force (g)');
hold off