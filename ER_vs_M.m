M = linspace(1, 4, 100);
ER = zeros(1, 100);
k = 1.257

for i = 1:100
   ER(i) = (1/M(i))*sqrt(((1+((k-1)/2)*M(i)^2)/(1+((k-1)/2)))^((k+1)/(k-1)));
end

ER_ideal = (1/3.08)*sqrt(((1+((k-1)/2)*3.08^2)/(1+((k-1)/2)))^((k+1)/(k-1)));

plot(M, ER)
hold on
plot(3.08, ER_ideal, 'rx');
xlabel('Mach number');
ylabel('Expansion Ratio');
hold off