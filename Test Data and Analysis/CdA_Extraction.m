%% Documentation
% Cd*A Extraction for Cold Flow Tests 1 & 2
% Original Author: Dan Zanko (09/04/2019)
clear, clc, close all

%% Setup
% loading test data files
load('ColdFlow1.mat')
load('ColdFlow2.mat')
% Constants
P_atm = 14.69594878; % assuming standard atmospheric pressure sea level (Laurel, MD is at approx. 150ft above sea level) [psi]
n_holes = 14; % number of orifices (holes) in showerhead injector plate tested
D_hole = 0.067; % diameter of each orifice hole [in]

% Startup Calcs
A_hole = pi*(D_hole^2)/4; % calculating discharge area of a single orifice [in^2]
A_inj = n_holes*A_hole/144; % total injector area (as designed) [in^2 --> ft^2]

%% Cold Flow 1
CF1.mdotsmthd = smooth(CF1.mdot);

CF1.dumptime = CF1.time(CF1.DumpInd(1,1):CF1.DumpInd(2,1)) - CF1.time(CF1.DumpInd(1,1)); % setting up a new time vector
CF1.dumpmdot = CF1.mdot(CF1.DumpInd(1):CF1.DumpInd(2));
CF1.dumpmdotsmthd = smooth(CF1.dumpmdot);

for n1 = CF1.DumpInd(1,1):CF1.DumpInd(2,1)
    ind = n1 - CF1.DumpInd(1) + 1;
    
    CF1.rho(ind,1) = N2Olookup("pressure",CF1.PT1(n1)*6.89476,0,'density')*0.00194032; % looks up density via PT1 and assuming Q = 0 and converting to US [slug/ft^3]
    CF1.Q(ind,1) = -(CF1.mdotsmthd(n1)*0.031081)/CF1.rho(ind,1);
    
    CF1.deltaP.Inj(ind,1) = (CF1.PT2(n1)); % gauge pressure (inj. - atm) [psig]
    CF1.CdA.Inj(ind,1) = abs((CF1.mdotsmthd(ind)*0.031081)/sqrt(abs(2*CF1.rho(ind,1)*CF1.deltaP.Inj(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF1.CdA.rawMdotInj(ind,1) = abs((CF1.mdot(ind)*0.031081)/sqrt(abs(2*CF1.rho(ind,1)*CF1.deltaP.Inj(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF1.Cd.Inj(ind,1) = CF1.CdA.Inj(ind,1)/A_inj;
    CF1.InjCavNum(ind,1) = (CF1.PT2(n1) - CF1.PT1(n1))/(CF1.PT2(n1)); % takes PT1 as sat pressure, discharge from PT2 to ambient

    CF1.deltaP.Tot(ind,1) = (CF1.PT1(n1)); % gauge pressure (tank - atm) [psig]
    CF1.CdA.Tot(ind,1) = abs((CF1.mdotsmthd(ind)*0.031081)/sqrt(abs(2*CF1.rho(ind,1)*CF1.deltaP.Tot(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF1.Cd.Tot(ind,1) = CF1.CdA.Tot(ind,1)/A_inj;
        
    CF1.deltaP.PTs(ind,1) = (CF1.PT1(n1)-CF1.PT2(n1)); % gauge pressure (tank - atm) [psig]
end

CF1.AVGmdot = mean(smooth(CF1.dumpmdotsmthd));

CF1.CdA.avgInj = mean(CF1.CdA.Inj(21:length(CF1.CdA.Inj)));
CF1.Cd.avgInj = mean(CF1.Cd.Inj(21:length(CF1.Cd.Inj)));

CF1.CdA.avgTot = mean(CF1.CdA.Tot(21:length(CF1.CdA.Tot)));
CF1.Cd.avgTot = mean(CF1.Cd.Tot(21:length(CF1.Cd.Tot)));

%% Cold Flow 2
CF2.mdotsmthd = smooth(CF2.mdot);

CF2.dumptime = CF2.time(CF2.DumpInd(1,1):CF2.DumpInd(2,1)) - CF2.time(CF2.DumpInd(1,1)); % setting up a new time vector
CF2.dumpmdot = CF2.mdot(CF2.DumpInd(1):CF2.DumpInd(2));
CF2.dumpmdotsmthd = smooth(CF2.dumpmdot);

for n2 = CF2.DumpInd(1,1):CF2.DumpInd(2,1)
    ind = n2 - CF2.DumpInd(1) + 1;
    
    CF2.rho(ind,1) = N2Olookup("pressure",CF2.PT2(n2)*6.89476,0,'density')*0.00194032; % looks up density via PT1 and assuming Q = 0 and converting to US [slug/ft^3]
    CF2.Q(ind,1) = -(CF2.mdotsmthd(ind)*0.031081)/CF2.rho(ind,1);
    
    CF2.deltaP.Inj(ind,1) = CF2.PT2(n2); % gauge pressure (inj. - atm) [psig]
    CF2.CdA.Inj(ind,1) = abs((CF2.mdotsmthd(ind)*0.031081)/sqrt(abs(2*CF2.rho(ind,1)*CF2.deltaP.Inj(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF2.Cd.Inj(ind,1) = CF2.CdA.Inj(ind,1)/A_inj;
    CF2.InjCavNum(ind,1) = (CF2.PT2(n2) - CF2.PT1(n2))/(CF2.PT2(n2)); % takes PT1 as sat pressure, discharge from PT2 to ambient

    CF2.deltaP.Tot(ind,1) = CF2.PT1(n2); % gauge pressure (tank - atm) [psig]
    CF2.CdA.Tot(ind,1) = abs((CF2.mdotsmthd(ind)*0.031081)/sqrt(abs(2*CF2.rho(ind,1)*CF2.deltaP.Tot(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF2.Cd.Tot(ind,1) = CF2.CdA.Tot(ind,1)/A_inj;
    
    CF2.deltaP.PTs(ind,1) = (CF2.PT1(n2)-CF2.PT2(n2)); % gauge pressure (tank - atm) [psig]
end

CF2.AVGmdot = mean(CF2.dumpmdotsmthd);

CF2.CdA.avgInj = mean(CF2.CdA.Inj(21:length(CF2.CdA.Inj)));
CF2.Cd.avgInj = mean(CF2.Cd.Inj(21:length(CF2.Cd.Inj)));

CF2.CdA.avgTot = mean(CF2.CdA.Tot(21:length(CF2.CdA.Tot)));
CF2.Cd.avgTot = mean(CF2.Cd.Tot(21:length(CF2.Cd.Tot)));

%% Pressure vs. Time Plots

figure, hold on
    plot(CF1.dumptime, CF1.PT1(CF1.DumpInd(1,1):CF1.DumpInd(2,1))),
    yyaxis left, plot(CF1.dumptime, CF1.PT2(CF1.DumpInd(1,1):CF1.DumpInd(2,1))),ylabel('Pressure (psig)')
    yyaxis right, plot(CF1.dumptime, CF1.InjCavNum),ylabel('Cavitation Number')
    xlabel('Time (s)'); grid on, grid minor
    title('Cold Flow 1 - Pressure and Cavitation Number vs. Time'), legend('Tank','Injector Manifold', 'Cavitation Number (injector)')
hold off

figure, hold on
    plot(CF2.dumptime, CF2.PT1(CF2.DumpInd(1,1):CF2.DumpInd(2,1))),
    yyaxis left, plot(CF2.dumptime, CF2.PT2(CF2.DumpInd(1,1):CF2.DumpInd(2,1))),ylabel('Pressure (psig)')
    yyaxis right, plot(CF2.dumptime, CF2.InjCavNum),ylabel('Cavitation Number')
    xlabel('Time (s)'); grid on, grid minor
    title('Cold Flow 2 - Pressure and Cavitation Number vs. Time'), legend('Tank','Injector Manifold', 'Cavitation Number (injector)')
hold off

figure, hold on
    subplot(2,1,1); hold on
        plot(CF1.dumptime, CF1.PT1(CF1.DumpInd(1,1):CF1.DumpInd(2,1)))
        plot(CF1.dumptime, CF1.PT2(CF1.DumpInd(1,1):CF1.DumpInd(2,1)))
        xlabel('Time (s)'), ylabel('Pressure (psig)'), grid on, grid minor
        title('Cold Flow 1 - Pressure vs. Time'),legend('Tank','Injector Manifold')
    hold off
    
    subplot(2,1,2); hold on
        plot(CF2.dumptime, CF2.PT1(CF2.DumpInd(1,1):CF2.DumpInd(2,1)))
        plot(CF2.dumptime, CF2.PT2(CF2.DumpInd(1,1):CF2.DumpInd(2,1)))
        xlabel('Time (s)'), ylabel('Pressure (psig)'),  grid on, grid minor
        title('Cold Flow 2 - Pressure vs. Time'),legend('Tank','Injector Manifold')
    hold off

hold off

%% mdot vs time plots
figure, hold on
    

%% CdA vs. DeltaP Plots
figure, hold on
    plot(CF1.deltaP.Tot(30:length(CF1.deltaP.Tot)), CF1.CdA.Inj(30:length(CF1.deltaP.Tot)))
    plot(CF2.deltaP.Tot(20:length(CF2.deltaP.Tot)), CF2.CdA.Inj(20:length(CF2.deltaP.Tot)))
    xlabel('\Delta P (psi)'), ylabel('C_DA (ft^2)'),  grid on, grid minor
    title('CdA vs. Total Pressure Drop'),legend('CF1','CF2')
hold off

%% mdot vs. CdA plots
figure, hold on
    plot(CF1.CdA.Inj(50:length(CF1.deltaP.Tot)),abs(CF1.mdotsmthd(50:length(CF1.deltaP.Tot))))
    plot(CF2.CdA.Inj(50:length(CF2.deltaP.Tot)),abs(CF2.mdotsmthd(50:length(CF2.deltaP.Tot))))
    xlabel('C_DA (ft^2)'), ylabel('Mass Flow Rate (lbm/s)'), grid on, grid minor
    title('Mass Flow Rate vs. C_DA'),legend('CF1','CF2')
hold off
    
%% Cleanup
misc.n1 = n1;
misc.n2 = n2;
misc.ind = ind;
clear n1 n2 ind