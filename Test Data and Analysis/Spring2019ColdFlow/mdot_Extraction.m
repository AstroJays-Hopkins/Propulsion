%% Documentation
% mdot Extraction for Cold Flow Tests 1 & 2
% Original Author: Dan Zanko (09/04/2019)

%% Cold Flow 1
clear
load('ColdFlow1.mat')
% getting mdot
for i = 2:length(CF1.mass)
    CF1.mdot(i) = (CF1.mass(i) - CF1.mass(i-1))/(CF1.time(i) - CF1.time(i-1));
end
clear i
save('ColdFlow1.mat')

%% Cold Flow 2
clear
load('ColdFlow2.mat')
for i = 2:length(CF2.mass)
    CF2.mdot(i) = (CF2.mass(i) - CF2.mass(i-1))/(CF2.time(i) - CF2.time(i-1));
end
clear i
save('ColdFlow2.mat')
clear