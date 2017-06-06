% Kuramoto-Sivashinsky equation
%
%         u_t = -u_xx - u_xxxx - (u^2/2)_x
%
% on [-16pi,16pi] with periodic boundary conditions.
% Long waves grow because of -u_xx;
% short waves decay because of -u_xxxx;
% the nonlinear term transfers energy from long to short.

% Grid, initial data, and plotting setup:
npts = 400;
h = 32*pi/npts;
x = -16*pi + (1:npts)'*h;
u = cos(x/16).*(1+sin(x/16));

Hf = figure(1); clf; hold on;
H = get(Hf, 'children');  set(H, 'fontsize', 16);
plt = plot(x, u, 'linewidth', 4);
axis([-16*pi 16*pi -4 4]);
grid on;
xlabel('x'); ylabel('u');


[I, Dxx, Dxc, Dxb, Dxf] = diff_matrices1d(length(x), h, 'p');

% First deriv, Laplacian and Biharmonic
D = Dxc;
L = Dxx;
H = L^2;


k = 0.1; %0.5*h;
Tf = 200;
numsteps = ceil(Tf/k);
%k = Tf/numsteps
Tf = k*numsteps;

plotgap = 2;
data = zeros(ceil(numsteps/plotgap) + 1, npts);
data(1, :) = u;

% matrix for implicit time-stepping
A = I + k*H + k*L;

%disp('type <return> to see solution'), pause

% Time-stepping:
t = 0;
for n=1:numsteps
  u = A\(u - k*(D*(u.^2/2)));
  t = n*k;
  %set(plt, 'ydata', u)
  %title(['t = ' num2str(t, 3)])
  %drawnow
  if mod(n, plotgap) == 0, data(n/plotgap + 1, :) = u; end
end

[maxs, ind1] = max(data, [], 2);
[mins, ind2] = min(data, [], 2);

fig2 = figure(2);
surf(x, 0:k*plotgap:Tf, data, 'Edgecolor', 'none'), hold on
plot3(ind1*h-16*pi, 0:k*plotgap:Tf, maxs, 'r.', 'markersize', 15)
%plot3(ind1*h-16*pi, 0:k*plotgap:Tf, mins, 'b.', 'markersize', 15)
hold off, axis([-16*pi 16*pi 0 Tf]), xlabel('Space'), ylabel('Time'),
title('Kuramoto-Sivashinsky'), colorbar, view([-10, 80]),
colormap(fig2, 'summer')

print('-bestfit', '../../Assignments/Assignment1/KS_plot', '-dpdf')

fprintf('Maximum is %1.4f, and occurs at x = %1.1f (%1.0f-th entry in the array).\n',...
    maxs(end), ind1(end)*h, ind1(end))

fprintf('Minimum is %1.4f, and occurs at x = %1.1f (%1.0f-th entry in the array).\n',...
    mins(end), ind2(end)*h, ind2(end))













