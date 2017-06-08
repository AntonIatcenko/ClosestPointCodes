% Question 1: demostratig the closest point computation for a circle.
R = 2;
numpoints = 50;
point = zeros(2, 2);
[Cx, Cy] = pol2cart(linspace(0, 2*pi, 64), R);
figure(1)
plot(Cx, Cy, 'linewidth', 2), hold on
axis([-R*3/2 R*3/2 -R*3/2 R*3/2]); axis square
title(['Closest points for the circle with radius ', num2str(R, 1)],...
    'Fontsize', 14)
for j = 1:numpoints
    point(1, :) = 3*R*rand([2 1]) - 3*R/2;
    point(2, :) = cp_circle(point(1, :), R);
    figure(1), plot(point(:, 1), point(:, 2), '-*', 'linewidth', 2)
end, hold off
function cpV = cp_circle(V, R), cpV = R/norm(V, 2)*V; end    