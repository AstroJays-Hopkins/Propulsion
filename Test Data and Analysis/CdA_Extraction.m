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
    CF1.deltaP(ind,1) = (CF1.PT1(n1) - P_atm); % gauge pressure [psig]
    CF1.CdA(ind,1) = abs((CF1.mdot(n1)*0.031081)/sqrt(abs(2*CF1.rho(ind,1)*CF1.deltaP(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
end
CF1.avgCdA = mean(CF1.CdA);
CF1.Cd = CF1.CdA/A_inj;
CF1.avgCd = mean(CF1.Cd);

%% Cold Flow 2
CF2.dumptime = CF2.time(CF2.DumpInd(1,1):CF2.DumpInd(2,1)) - CF2.time(CF2.DumpInd(1,1)); % setting up a new time vector

for n2 = CF2.DumpInd(1,1):CF2.DumpInd(2,1)
    ind = n2 - CF2.DumpInd(1) + 1;
    CF2.rho(ind,1) = N2Olookup("pressure",CF2.PT2(n2)*6.89476,0,'density')*0.00194032; % looks up density via PT1 and assuming Q = 0 and converting to US [sl;ug/ft^3]
    CF2.deltaP(ind,1) = (CF2.PT2(n2) - P_atm); % gauge pressure [psig]
    CF2.CdA(ind,1) = abs((CF2.mdot(n2)*0.031081)/sqrt(abs(2*CF2.rho(ind,1)*CF2.deltaP(ind,1)*144))); % calculating Cd*A (mdot is converted from lbm/s to slug/s) (deltaP is converted from psig to psfg)
end
CF2.avgCdA = mean(CF2.CdA);
CF2.Cd = CF2.CdA/A_inj;
CF2.avgCd = mean(CF2.Cd);