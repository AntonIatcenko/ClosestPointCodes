%% Radial Basis Functions Method from Scratch
%
%% Grid
N = 10;   Ne = 50;        % Sizes of sampling and plotting grid
x = rand(N, 1);           % Sampling grid
xe = linspace(0, 1, Ne);  % Plotting grid
fun = @(x) sin(2*x);      % Underlying function
f = fun(x);               % Grid values
%% Radial Basis Interpolation
phi = @(ep, r) exp(-ep^2*r.^2);
eps = 1;

[X1, X2] = meshgrid(x);

%%
A = phi(eps, abs(x-x'));
%A = phi(eps, abs(X2 - X1));
lam = A\f;

%%
B = phi(eps, x - xe)';

fint = B*lam;

%% Plot                       
figure(1)
plot(x, f, 'r.', 'markersize', 30), hold on
plot(xe, fint, 'linewidth', 3)
plot(xe, fun(xe), '-', 'linewidth', 1), hold off
title(['RBF interpolation with Gaussian functions and eps = ',...
    num2str(eps, 1)], 'fontsize', 16)
legend({'Sampled data', 'Interpolant', 'Underlying function'},...
    'location', 'northwest', 'fontsize', 14)