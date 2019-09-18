%% Documentation
% Cd*A Extraction for Cold Flow Tests 1 & 2
% Original Author: Dan Zanko (09/04/2019)
clear, clc
%% Setup
% loading test data files
load('ColdFlow1.mat')
load('ColdFlow2.mat')
% Constants
P_atm = 14.69594878; % assuming standard atmospheric pressure sea level (Laurel, MD is at approx. 150ft above sea level) [psi]
A_inj = (0.07221642243)/144; % injector area (as designed) [in^2 --> ft^2]

%% Cold Flow 1
CF1.dumptime = CF1.time(CF1.DumpInd(1,1):CF1.DumpInd(2,1)) - CF1.time(CF1.DumpInd(1,1)); % setting up a new time vector

for n1 = CF1.DumpInd(1,1):CF1.DumpInd(2,1)
    ind = n1 - CF1.DumpInd(1) + 1;
    CF1.rho(ind,1) = N2Olookup("pressure",CF1.PT1(n1)*6.89476,0,'density')*0.00194032; % looks up density via PT1 and assuming Q = 0 and converting to US [slug/ft^3]
    CF1.Q(ind,1) = -(CF1.mdot(n1)*0.031081)/CF1.rho(ind,1);
    CF1.deltaP.Inj(ind,1) = (CF1.PT2(n1)); % gauge pressure (inj. - atm) [psig]
    CF1.CdA.Inj(ind,1) = abs((CF1.mdot(n1)*0.031081)/sqrt(abs(2*CF1.rho(ind,1)*CF1.deltaP.Inj(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF1.Cd.Inj(ind,1) = CF1.CdA.Inj(ind,1)/A_inj;
    CF1.deltaP.Tot(ind,1) = (CF1.PT1(n1)); % gauge pressure (tank - atm) [psig]
    CF1.CdA.Tot(ind,1) = abs((CF1.mdot(n1)*0.031081)/sqrt(abs(2*CF1.rho(ind,1)*CF1.deltaP.Tot(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF1.Cd.Tot(ind,1) = CF1.CdA.Tot(ind,1)/A_inj;
    CF1.InjCavNum(ind,1) = (CF1.PT2(n1) - CF1.PT1(n1))/(CF1.PT2(n1)); % takes PT1 as sat pressure, discharge from PT2 to ambient
    CF1.deltaP.PTs(ind,1) = (CF1.PT1(n1)-CF1.PT2(n1)); % gauge pressure (tank - atm) [psig]
   

end

CF1.CdA.avgInj = mean(CF1.CdA.Inj(21:length(CF1.CdA.Inj)));
CF1.Cd.avgInj = mean(CF1.Cd.Inj(21:length(CF1.Cd.Inj)));
CF1.CdA.avgTot = mean(CF1.CdA.Tot(21:length(CF1.CdA.Tot)));
CF1.Cd.avgTot = mean(CF1.Cd.Tot(21:length(CF1.Cd.Tot)));

figure, hold on
plot(CF1.dumptime, CF1.PT1(CF1.DumpInd(1,1):CF1.DumpInd(2,1))),
[CF1Ax,CF1Line1,CF1Line2] = plotyy(CF1.dumptime, CF1.PT2(CF1.DumpInd(1,1):CF1.DumpInd(2,1)), CF1.dumptime, CF1.InjCavNum)
ylabel(CF1Ax(1),'Pressure (psig)')
ylabel(CF1Ax(2),'Cavitation Number')

CF1Line2.LineStyle = ':';
legend('PT1','PT2')


%% Cold Flow 2
CF2.dumptime = CF2.time(CF2.DumpInd(1,1):CF2.DumpInd(2,1)) - CF2.time(CF2.DumpInd(1,1)); % setting up a new time vector

for n2 = CF2.DumpInd(1,1):CF2.DumpInd(2,1)
    ind = n2 - CF2.DumpInd(1) + 1;
    CF2.rho(ind,1) = N2Olookup("pressure",CF2.PT2(n2)*6.89476,0,'density')*0.00194032; % looks up density via PT1 and assuming Q = 0 and converting to US [slug/ft^3]
    CF2.Q(ind,1) = -(CF2.mdot(n2)*0.031081)/CF2.rho(ind,1);
    CF2.deltaP.Inj(ind,1) = CF2.PT2(n2); % gauge pressure (inj. - atm) [psig]
    CF2.CdA.Inj(ind,1) = abs((CF2.mdot(n2)*0.031081)/sqrt(abs(2*CF2.rho(ind,1)*CF2.deltaP.Inj(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF2.Cd.Inj(ind,1) = CF2.CdA.Inj(ind,1)/A_inj;
    CF2.deltaP.Tot(ind,1) = CF2.PT1(n2); % gauge pressure (tank - atm) [psig]
    CF2.CdA.Tot(ind,1) = abs((CF2.mdot(n2)*0.031081)/sqrt(abs(2*CF2.rho(ind,1)*CF2.deltaP.Tot(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
    CF2.Cd.Tot(ind,1) = CF2.CdA.Tot(ind,1)/A_inj;
end
CF2.CdA.avgInj = mean(CF2.CdA.Inj(21:length(CF2.CdA.Inj)));
CF2.Cd.avgInj = mean(CF2.Cd.Inj(21:length(CF2.Cd.Inj)));
CF2.CdA.avgTot = mean(CF2.CdA.Tot(21:length(CF2.CdA.Tot)));
CF2.Cd.avgTot = mean(CF2.Cd.Tot(21:length(CF2.Cd.Tot)));

