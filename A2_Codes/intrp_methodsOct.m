%% Spatial Convergence Study      
% for Laplace-Beltramy on a circle of radius R:
% Lap*u = f
%% Parameters                     
R = 1;    N = 10;       % Radius and number of trials
grids = 2.^(3:N+2);     % Grids
errors = zeros(4, N);   % Preallocation
m = {'nearest', 'linear', 'cubic', 'pchip'};
%% Iterations   
tic
for q=1:4
for j = 1:N
    dx = 3*R/(grids(j));              % Resolution
    x = dx-3*R/2:dx:3*R/2;          % 1d grid
    [X, Y] = meshgrid(x, x);        % 2d grid
    [th, ~] = cart2pol(X, Y);       % Circle
    [cpx, cpy] = pol2cart(th, R);   % Closest points
    u = sin(4*th);                  % Function to be differentiated
    I = speye(grids(j));
    D = (circshift(I, [1, 0]) - 2*I + circshift(I, [-1, 0]))/(dx^2);  % 2nd derivative matrix
    LapU = D*u + u*D;               % Derivative on the grid
    f = interp2(X, Y, LapU, cpx, cpy, m{q});
    errors(q, j) = norm(f(:) + 4^2*u(:), inf);
    %fprintf('N = %03i, error = %2.4f \n', grids(j), supError)
end
end
%% Convergence Rate   
p1 = polyfit(log(grids), log(errors(1, :)), 1);  
p2 = polyfit(log(grids), log(errors(2, :)), 1);
p3 = polyfit(log(grids), log(errors(3, :)), 1);
p4 = polyfit(log(grids), log(errors(4, :)), 1);
pp = [p1 p2 p3 p4];   rates = -[p1(1) p2(1) p3(1) p4(1)];
%% Plotting 
figure(2)
loglog(grids, errors(1, :), 'r-s', 'linewidth', 2, 'markersize', 10), hold on
loglog(grids, errors(2, :), 'b-s', 'linewidth', 2, 'markersize', 10)
loglog(grids, errors(3, :), 'c-s', 'linewidth', 2, 'markersize', 10)
loglog(grids, errors(4, :), 'k-s', 'linewidth', 2, 'markersize', 10), hold off
title('Spatial convergence rates for different interpolation methods', 'fontsize', 16)
legend(['nearest - ', num2str(rates(1))],...
          ['linear - ' , num2str(rates(2))],...
          ['cubic - '  , num2str(rates(3))],...
          ['pchip - '  , num2str(rates(4))])                        

    

    

