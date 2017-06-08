%% Spatial Convergence Study      
% for Laplace-Beltramy on a circle of radius R:
% Lap*u = f
%% Parameters                     
R = 1;    N = 10;       % Radius and number of trials
grids = 2.^(3:N+2);     % Grids
errors = zeros(1, N);   % Preallocation
%% Iterations                     
for j = 1:N
    dx = 3*R/grids(j);              % Resolution
    x = dx-3*R/2:dx:3*R/2;          % 1d grid
    [X, Y] = meshgrid(x, x);        % 2d grid
    [th, ~] = cart2pol(X, Y);       % Circle
    [cpx, cpy] = pol2cart(th, R);   % Closest points
    u = sin(2*th);                  % Function to be differentiated
    I = speye(grids(j));
    D = (circshift(I, [1, 0]) - 2*I + circshift(I, [-1, 0]))/(dx^2);  % 2nd derivative matrix
    LapU = D*u + u*D;               % Derivative on the grid
    f = interp2(X, Y, LapU, cpx, cpy, 'cubic');
    errors(j) = norm(f(:) + 2^2*u(:), inf);
    %fprintf('N = %03i, error = %2.4f \n', grids(j), supError)
end
%% Convergence Rate               
p = polyfit(log(grids), log(errors), 1);  % Fitting a linear model into 
rate = -p(1);                             % logarithmically scaled data
%% Plotting                       
figure(1)
loglog(grids, exp( polyval(p, log(grids))), 'linewidth', 2);
hold on, loglog(grids, errors, '.', 'markersize', 20); hold off
title(['Spatial convergence rate is ', num2str(rate, 3)], 'fontsize', 16)
legend({'Linear model', 'Errors'}, 'fontsize', 16)
print('-bestfit', '../../Assignments/Assignment2/HeatCircleConv1', '-dpdf')
%%
%% Temporal Convergence Study     
% for heat equation on a circle of radius R
%% Parameters                     
R = 1;    N = 9;           % Radius and number of trials
M = 100;  dx = 3*R/M;      % Spatial parameters
x = dx-3*R/2:dx:3*R/2;
Tfinal = 1;                % Length of the simulations
Tsteps = 2.^(5:N+4);       % Temporal resolutions
errors = zeros(1, N);      % Preallocation
m = 2;                     % Fourier mode number for IC
nu = 1;                    % Diffusivity
[X, Y] = meshgrid(x, x);       
[th, ~] = cart2pol(X, Y);   
[cpx, cpy] = pol2cart(th, R); 
u0 = sin(m*th(:));  
u_true = exp(-Tfinal*nu*m^2)*u0; 
I = speye(M);  II = speye(M^2);
D = (circshift(I, [1, 0]) - 2*I + circshift(I, [-1, 0]))/(dx^2);        
Lap = nu*(kron(I, D) + kron(D, I));   
%% Iterations                     
parfor j = 1:N
    dt = Tfinal/Tsteps(j);  [Low, Up] = lu(II - dt*Lap); u = u0;
    for t = 1:Tsteps(j)  
        u = Up\(Low\u);             % Time step on the embedding grid   
        u = interp2(X, Y, reshape(u, M, M), cpx, cpy, 'cubic');   % Interpolation to the circle
        u = u(:); 
    end
    errors(j) = norm(u - u_true, inf);
end
%% Convergence Rate               
p = polyfit(log(Tsteps), log(errors), 1);  % Fitting a linear model into 
rate = -p(1);                             % logarithmically scaled data
%% Plotting                       
figure(2)
loglog(Tsteps, exp( polyval(p, log(Tsteps))), 'linewidth', 2);
hold on, loglog(Tsteps, errors, '.', 'markersize', 20); hold off
title(['Temporal convergence rate is ', num2str(rate, 3)], 'fontsize', 16)
legend({'Linear model', 'Errors'}, 'fontsize', 16)                        
print('-bestfit', '../../Assignments/Assignment2/HeatCircleConv2', '-dpdf')