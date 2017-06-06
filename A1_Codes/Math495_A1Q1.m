%% Convergence Study for the Second Order Finite Difference Sheme
%% Set up                           
M = 20;                 % Number of trials
errors = zeros(1, M);   % Preallocating for errors
%% Running Trials                   
for G = 1:M
    N = 2^G;                  % Grid size
    dx = 1/N;                 % Spatial resolution
    x = 0:dx:1;               % Grid
    f = sin(pi*x/2);          % Function sampled at the grid
    fx = pi/2*cos(pi*x/2);    % Analytical derivative 

    % Derivative on the interior using centered difference 
    fxnum = ([f(2:end) f(1)] - [f(end) f(1:end-1)])/dx/2; 
    % Derivative at x = 0 using forward difference
    fxnum(1) = -(3*f(1) - 4*f(2) + f(3))/dx/2;
    % Derivative at x = 1 using backward difference
    fxnum(end) = (3*f(end) - 4*f(end-1) + f(end-2))/dx/2;
    % Error
    errors(G) = sqrt(dx)*norm( fx-fxnum, 2);
end
%% Convergence Analysis             
ConvRate = polyfit(log(2)*(1:M-2), log(errors(1:end-2)), 1);
figure(1)
loglog(2.^(1:M), exp(polyval(ConvRate, log(2)*(1:M))), 'linewidth', 3), hold on
loglog(2.^(1:M), errors, '.', 'markersize', 20), hold off
title(['Convergence Rate is ', num2str( -ConvRate(1), 3)], 'fontsize', 14)
legend({'Errors', 'Best Fit Linear Model'}, 'Fontsize', 14)
xlabel('log(N)'), ylabel('log( E_2 )')
%% Outputs                          
print('-bestfit', '../../Assignments/Assignment1/ConvRate', '-dpdf')
table = [1./(2.^(1:G) - 1)' errors'];
dlmwrite('../../Assignments/Assignment1/table.csv', table, 'precision', '%2.2g')