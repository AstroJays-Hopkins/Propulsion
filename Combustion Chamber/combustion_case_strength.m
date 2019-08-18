P = 3.6e6; %[Pa]
SF_min = 3;
r_in = 3.125*0.0254/2; %[m] 
r_out = 3.5*0.0254/2; %[m]

r = (r_out+r_in)/2

t = 0.1875*0.0254;

T = [24, 100, 149, 204, 260];
y_strength = [310e6, 290e6, 234e6, 131e6, 51e6]; 

stress = P * r / t;

SF = y_strength/stress;

plot(T, SF)
hold on
xlabel('Temp (C)')
ylabel('SF')
plot(T, SF_min*ones(1, length(T)))
hold off