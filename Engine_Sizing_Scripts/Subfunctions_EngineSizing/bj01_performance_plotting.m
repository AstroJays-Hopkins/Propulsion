% Plotting for the BJ01 Performance Analysis Script

clear

load BJ01_motor_params_102218.mat
% 
% figure(1)
% hold on
% plot((1:t_burn/deltat)*deltat, P_RT/PSItoPa, 'k')
% axis([0, 13, 0, 850]);
% title('Run Tank Pressure');
% xlabel('Time (s)')
% ylabel('Pressure (psi)');
% hold off
% 
% figure(2)
% plot((1:t_burn/deltat)*deltat, mdot_ox*KGtoLBM, 'k')
% hold on
% axis([0, 13, 0, 4]);
% title('N2O Mass Flow Rate');
% xlabel('Time (s)')
% ylabel('Mass Flow Rate (lbm/s)');
% hold off
% 
% figure(3)
% plot((1:t_burn/deltat)*deltat, P_c/PSItoPa, 'k')
% hold on
% axis([0, 13, 0, 1000]);
% title('Combustion Chamber Pressure');
% xlabel('Time (s)')
% ylabel('Pressure (psi)');
% hold off



% figure(4)
% plot((1:t_burn/deltat)*deltat, OFR, 'k');
% axis([0, 13, 0, 12]);
% title('OF Ratio');
% xlabel('Time (s)');
% ylabel('Oxidizer to Fuel Ratio')

% for t = 1:((t_burn/deltat)-1)
%     [T_c_t(t), k_c_t, R_t] = combustionParameters(OFR(t), P_c(t));
% end
% T_c_t(t+1) = T_c_t(t);
% 
% 
% figure(5)
% plot((1:t_burn/deltat)*deltat, T_c_t, 'k');
% hold on
% axis([0, 13, 3100, 3350]);
% title('Combustion Temperature');
% xlabel('Time (s)');
% ylabel('Temperature (K)')
% hold off
% 
figure(6)
plot((1:t_burn/deltat)*deltat, rdot*MtoIN, 'k');
hold on
axis([0, 13, 0, 0.07]);
title('Fuel Regression Rate')
xlabel('Time (s)');
ylabel('Regression rate (in/s)')
hold off

% 
% figure(7)
% hold on
% plot((1:t_burn/deltat)*deltat, T/LBFtoN, 'k');
% axis([0, 13, 0, 1000]);
% title('BJ-01 Thrust vs Time');
% xlabel('Time (s)');
% ylabel('Thrust (lbf)')