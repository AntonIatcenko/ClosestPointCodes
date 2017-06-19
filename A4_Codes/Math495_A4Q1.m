%% Interpolation with Radial Basis Functions
%
%
%% Parameters
Nsample = 41;         % Number of samples
Nint = 5e2;           % Size of interpolation grid
%% Grids
H = halton(41, 2);      % Sampling grid
dx = 1/(Nint-1);        % Spacing for thse interpolation grid
x = 0:dx:1;             % 1d interpolation grid
[X, Y] = meshgrid(x);   % 2d interpolation grid
%% Data
fk = 59./(67 + (H(:, 1) + 1/7).^2 + (H(:, 2) - 1/11).^2);  % Data to be interpolated
ftrue = 59./(67 + (X(:) + 1/7).^2 + (Y(:) - 1/11).^2);     % True values 
%% Radial Basis Interpolation
phi = @(ep, r) exp(-ep^2*r);
dist1 = ( H(:, 1) - H(:, 1)' ).^2 + ( H(:, 2) - H(:, 2)' ).^2;
dist2 = ( H(:, 1)' - X(:) ).^2 + ( H(:, 2)' - Y(:) ).^2;
%% Iterations
eps_vals = 1; %linspace(0.05, 1, 32);
L = length(eps_vals);
errors = zeros(1, L);
condnums = zeros(1, L);

warning('off', 'MATLAB:nearlySingularMatrix')
tic
for j = 1:L

    eps = eps_vals(j);
    A = phi(eps, dist1);
    B = phi(eps, dist2);
    lam = A\fk;
    fapprox = B*lam;
    errors(j) = norm(fapprox - ftrue, inf);
    condnums(j) = cond(A);
    
end
toc
warning('on', 'MATLAB:nearlySingularMatrix')
%%

%%
figure(1)
semilogy(eps_vals, errors, '.', 'markersize', 20)
xlabel('\epsilon', 'fontsize', 20)
ylabel('Error', 'fontsize', 20)
title('Errors', 'fontsize', 14)

figure(2)
semilogy(eps_vals, condnums, '.', 'markersize', 20)
title('Condition numbers', 'fontsize', 14)

%
%%
% fapproxvec = reshape(fapprox, Nint, Nint);
% 
% figure(1)
% surf(X, Y, fapproxvec, 'edgecolor', 'none'), hold on
% plot3(H(:, 1), H(:, 2), fk, 'r.', 'markersize', 20)
% 
% fprintf('Maximum error is %2.4f.\n', norm(fapproxvec - ftrue, inf))
% 






























