%% Advection equation on a circle     
% convergence study
%% Parameters                         
M = 05;    Tfinal = 1;          % Number of trials
grids = 2.^(6:M+5);             % Spatial grids of different sizes
errors = zeros(1, M);           % Preallocating for errors
uex = @(th,t) (cos(th-t)).^3;   % Analytic solution
%
%% Trial Runs                         
parfor j = 1:M
dx = 1/grids(j);     x1d = (-2.0:dx:2.0)';    [x, y] = meshgrid(x1d, x1d);
[cpx, cpy, dist] = cpCircle(x, y);   dim = 2;  p = 3;  
bw = rm_bandwidth(dim, p);   band = find(abs(dist) <= bw*dx);
cpx = cpx(band); cpy = cpy(band);     x = x(band); y = y(band);
[th, ~] = cart2pol(cpx, cpy);   u = uex(th, 0);  
w1 = -sin(th);   w2 = cos(th);
E = interp2_matrix(x1d, x1d, cpx, cpy, p, band);
L = laplacian_2d_matrix(x1d, x1d, 2, band);
[Dxb, Dxf, Dyb, Dyf] = firstderiv_upw1_2d_matrices(x1d, x1d, band);
dt = 0.25*dx;  numsteps = ceil(Tfinal/dt);  dt = Tfinal / numsteps;
    for kt = 1:numsteps
        rhs = -( (w1 < 0) .* (Dxf*(u.*w1)) + (w1 >= 0) .* (Dxb*(u.*w1)) + ...
           (w2 < 0) .* (Dyf*(u.*w2)) + (w2 >= 0) .* (Dyb*(u.*w2)) );
          unew = u + dt*rhs;  u = E*unew;
    end
    errors(j) = norm(u - uex(th, Tfinal), inf);
end
%% Convergence Analysis               
p = polyfit(log(grids), log(errors), 1);
figure(1)
loglog(grids, exp( polyval(p, log(grids))), 'linewidth', 2) 
hold on
loglog(grids, errors, '.', 'markersize', 20)
title(sprintf('Convergence rate is %1.2f', -p(1)), 'fontsize', 20)
legend({'Best fit linear model', 'Errors'}, 'fontsize', 16)