hole_pos = zeros(15, 2);
r = 2.5; %[in] bolt diameter radius
spacing = 360/15; %[deg]

for i = 1:(length(hole_pos))
    hole_pos(i, 1) = r*sind(spacing*(i-1));
    hole_pos(i, 2) = r*cosd(spacing*(i-1));
end

plot(hole_pos(:,1), hole_pos(:,2), 'x')