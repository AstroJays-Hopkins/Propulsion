%unit conversions
PSItoPa = 6894.74; %[psi to Pa]
LBFtoN = 4.44822; %[lbf to N]
MtoIN = 39.3701; %[m to in]
KGtoLBM = 2.20462; %[kg to lbm]

load('BJ01_motor_params_040819.mat');

fid = fopen('BJ-01_motor_params_040819.txt', 'w');

fprintf(fid, 'ASTROJAYS\nBJ-01 MOTOR PARAMETERS CALCULATED BY thrust_curve.m\nUpdated Mar 16, 2019\n\n\n\n');

fprintf(fid, 'GENERAL MOTOR PERFORMANCE\n');
fprintf(fid,'**************************************************************\n');
fprintf(fid, 'Isp (s): %f\tTotal Impulse (N.s); %f\nBurn Time (s): %f\tMax Thrust (lbf): %f\nOF Ratio: %f\nFlight Ceiling, VERY rough esimate (ft): %f\nMaximum G Force (g): %f\n', Isp, I_total, t_burn, max(T(3:end))/LBFtoN, OF, max(altitude_FS)*MtoIN/12, max(g_force(3:end)));
fprintf(fid,'**************************************************************\n\n\n');

fprintf(fid, 'INITIAL ISENTROPIC FLOW CALCULATIONS\n');
fprintf(fid, '**************************************************************\n');
fprintf(fid, 'CHAMBER CONDITIONS:\n\n');
fprintf(fid, 'Pressure (psi): %f\tTemperature (K): %f\nVelocity (m/s): 0\n**************************************************************\n', P_c(20)/PSItoPa, T_c);
fprintf(fid, 'THROAT CONDITIONS:\n\n');
fprintf(fid, 'Pressure (psi): %f\tTemperature (K): %f\nVelocity (m/s): %f\tDensity (kg/m^3): %f\n**************************************************************\n', P_star/PSItoPa, T_star, u_star, rho_star);
fprintf(fid, 'NOZZLE EXIT CONDITIONS:\n\n');
fprintf(fid, 'Pressure (psi): %f\tTemperature (K): %f\nVelocity (m/s): %f\n**************************************************************\n\n\n\n', P_e/PSItoPa, T_e, u_e(20));

fprintf(fid, 'ENGINE DESIGN PARAMETERS\n');
fprintf(fid, '**************************************************************\n');
fprintf(fid, 'Propellant Mass (lbm): %f\n', m_fuel*KGtoLBM);
fprintf(fid, 'Initial Inner Port Radius (in): %f\nOuter Propellant Radius (in): %f\nLength of Propellant Grain (in): %f\n\n', rin_fuel*MtoIN, rout_fuel*MtoIN, L_fuel*MtoIN);
fprintf(fid, 'Expansion Ratio: %f\nInlet Diameter (in): %f\nThroat Diameter (in): %f\nExit Diameter (in): %f\n', ER, rout_fuel*2*MtoIN, D_star*MtoIN, D_e*MtoIN);
fprintf(fid, '**************************************************************\n\n\n\n');

fprintf(fid, 'RUN TANK DESIGN PARAMETERS\n');
fprintf(fid, '**************************************************************\n');
fprintf(fid, 'Oxizider Mass (lbm): %f\nEstimated Ox Mass Flow Rate (lbm/s): %f\n', m_ox*KGtoLBM, mdot_ox(20)*KGtoLBM);
fprintf(fid, '**************************************************************\n\n\n\nEND OUTPUT');

fclose(fid);