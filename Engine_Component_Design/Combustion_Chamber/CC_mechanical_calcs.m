%%SF calculator for CC PV
% Written by Andrew Colombo
% Astrojays Summer 2018
% Last edit: 8 Aug 2018


%% Geometry
%Nozzle
r_i = 1.5;
r_t = 0.515;
r_e = 1.256;

A_i = pi*(r_i^2 - r_t^2); %Area, inlet
A_e = pi*(r_e^2 - r_t^2); %Area, exit

%CC Case
r_cc_o = 2.25; %[in] Outer radius of combustion chamber PV
r_cc_i = 2;

%% Loading
P_nom = 522; %[psi]
P_meop = P_nom * 1.5; %Max Expected Op Pressure
P_atm = 14; %[psi]

F_P = (P_meop - P_atm) * A_i; %[lbf] Force due to pressure

%% Mat Props
S_u_6061T6 = 45000; %[psi] Ultimate strength of 6061-T6
nu_6061T6 = 0.3;

%% Bolts Calcs, 1/4-20 High Strength Steel Hex Screws (https://www.mcmaster.com/#90002a104/=1e2r979)
n_bolts = 15;
F_bolt = F_P/n_bolts;

A_t = 0.031; %[in^2]
S_u_HSS = 150000; %[psi] Ultimate tensile strength of "High Strength Steel" as per McMaster-Carr
FS_bolts = S_u_HSS/(F_bolt/A_t);

%% PV Longitudinal Stress
FS_long = S_u_6061T6/(F_P/(pi*(r_cc_o^2 - r_cc_i^2)));

%% PV Hoop
FS_hoop_case = S_u_6061T6/(P_meop*((r_cc_i+r_cc_o)/2)/(r_cc_o-r_cc_i));

%% PV Hydrotest Cap - max stress derived from circular plate fixed at edges
FS_cap = 3;
t_cap = sqrt(FS_cap*3*P_meop*r_cc_o^2/(4*S_u_6061T6));

%% Flange, assuming weld material 
% HAZ and knockdown info: https://www.esabna.com/us/en/education/blog/the-haz-in-aluminum-welds.cfm
S_u_6061T6_HAZ = 27000; %[psi] strength of 6061T6 HAZ
S_shear_allowable = 0.3*27000; %Shigley's Table 9.4
r_bolts = 2.75; %[in] Radius of bolt hole centerlines
height_weld = 0.125; %[in]
A_weld_throat_area = 1.414 * height_weld * (2*pi*r_cc_o);
t_flange = 0.125;
I_butt = ((2*pi*r_cc_o)*height_weld^3)/12 + (t_flange/2)^2 * height_weld * (2*pi*r_cc_o);

M_weld = F_P * (r_bolts - r_cc_o);
tensile_stress_max_weld = M_weld*(t_flange/2)/I_butt;
shear_stress_max_weld = F_P/(height_weld * (2*pi*r_cc_o));
FS_weld_shear = S_u_6061T6_HAZ/shear_stress_max_weld;
FS_weld_tensile = S_u_6061T6_HAZ/tensile_stress_max_weld;


%% Nozzle Compression
S_u_comp_gr = 8100; %[psi]
A_crit_noz = pi*(1.15^2 - 1^2)
FS_gr_comp = S_u_comp_gr/(F_P/A_crit_noz)

%% Aft closeout, Table 11.2, entry 1e, Roark's Formulas for Stress and Strain, 7th ed
a_ac = 2.75;
r0_ac = 1.75;
b_ac = 1.55;
w_ac = F_P/(2*pi*r0_ac);
t_ac = 0.25;

L9 = (r0_ac/a_ac)*((1+nu_6061T6)/2*log(a_ac/r0_ac)+(1-nu_6061T6)*(1-(r0_ac/a_ac)^2)/4);
L6 = r0_ac*((r0_ac/a_ac)^2-1+2*log(a_ac/r0_ac))/(4*a_ac)
C4 = 0.5*((1+nu_6061T6)*(b_ac/a_ac)+(1-nu_6061T6)*(a_ac/b_ac))
C7 = 0.5*(1-nu_6061T6^2)*(a_ac/b_ac-b_ac/a_ac)

M_ac = w_ac*a_ac*(L9-C7*L6/C4);
tensile_stress_ac = M_ac*(t_ac/2)/(2*pi*b_ac*t_ac^3/12);

FS_closeout_tensile = S_u_6061T6/tensile_stress_ac;

%% Print FSs
fprintf('\n\nFactor of Safety:\n')
fprintf('Bolts:\t\t\t%f\n', FS_bolts)
fprintf('CC Hoop:\t\t%f\n', FS_hoop_case)
fprintf('CC Longitudinal:\t%f\n', FS_long)
fprintf('Flange Weld Shear:\t%f\n', FS_weld_shear)
fprintf('Flange Weld Tensile:\t%f\n', FS_weld_tensile)
fprintf('PV Cap:\t\t\t%f\n', FS_cap)
fprintf('Nozzle Compression:\t%f\n', FS_gr_comp)
fprintf('Aft Closeout:\t\t%f\n', FS_closeout_tensile)