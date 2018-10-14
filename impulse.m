im = 0;

for i = 1:1300
    im = im + F_T(i)*deltat;
end

im = im * 0.224809