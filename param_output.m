load('motor_param_output.mat');

fid = fopen('BJ-01_motor_params.txt', 'w');

fprintf(fid, 'ASTROJAYS\nBJ-01 MOTOR PARAMETERS CALCULATED BY thrust_curve.m\nUpdated Feb 22, 2018\n\n\n\n');

fprintf(fid, 'GENERAL MOTOR PERFORMANCE\n');
fprintf(fid,'**************************************************************\n');
fprintf(fid, 'Isp (s): %f\tTotal Impulse (N.s); %f\nBurn Time (s): %f\tMax Thrust (N): %f\nOF Ratio: %f\nFlight Ceiling, VERY rough esimate (m): %f\nMaximum G Force (g): %f\n', Isp, I_total, t_burn, Thrust_max, OF, max(altitude_FS), max(g_force));
fprintf(fid,'**************************************************************\n\n\n');

fprintf(fid, 'ISENTROPIC FLOW CALCULATIONS\n');
fprintf(fid, '**************************************************************\n');
fprintf(fid, 'CHAMBER CONDITIONS:\n\n');
fprintf(fid, 'Pressure (Pa): %f\tTemperature (K): %f\nVelocity (m/s): 0\n**************************************************************\n', P_c, T_c);
fprintf(fid, 'THROAT CONDITIONS:\n\n');
fprintf(fid, 'Pressure (Pa): %f\tTemperature (K): %f\nVelocity (m/s): %f\tDensity (kg/m^3): %f\n**************************************************************\n', P_star, T_star, u_star, rho_star);
fprintf(fid, 'NOZZLE EXIT CONDITIONS:\n\n');
fprintf(fid, 'Pressure (Pa): %f\tTemperature (K): %f\nVelocity (m/s): %f\n**************************************************************\n\n\n\n', P_e, T_e, u_e);

fprintf(fid, 'ENGINE DESIGN PARAMETERS\n');
fprintf(fid, '**************************************************************\n');
fprintf(fid, 'Propellant Mass (kg): %f\nOxizider Mass (kg): %f\nEstimated Ox Mass Flow Rate (kg/s): %f\n', m_fuel, m_ox, mdot_ox);
fprintf(fid, 'Initial Inner Port Radius (m): %f\nOuter Propellant Radius (m): %f\nLength of Propellant Grain (m): %f\n\n', rin_fuel, rout_fuel, L_fuel);
fprintf(fid, 'Expansion Ratio: %f\nInlet Diameter (m): %f\nThroat Diameter (m): %f\nExit Diameter (m): %f\n', ER, sqrt(4*D_star^2), D_star, D_e);
fprintf(fid, '**************************************************************\n\n\n\nEND OUTPUT');
fclose(fid);