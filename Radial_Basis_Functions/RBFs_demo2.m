%% Radial Basis Functions Method from Scratch
%
%% Grid
N = 100;   Ne = 100;                    % Sizes of sampling and plotting grid
x = rand(N, 2);                        % Sampling grid
temp = linspace(0, 1, Ne);             % Plotting grid
[xe, ye] = meshgrid(temp);
xe = xe(:);  ye = ye(:);
fun = @(x, y) sin(3*x).*cos(5*y);      % Underlying function
f = fun(x(:, 1), x(:, 2));             % Grid values
%% Radial Basis Interpolation
phi = @(ep, r) exp(-ep^2*r.^2);
eps = 1;

%%
[X1, X2] = meshgrid(x(:, 1));
[Y1, Y2] = meshgrid(x(:, 2));
[Xe1, Xe2] = meshgrid( x(:, 1), xe);
[Ye1, Ye2] = meshgrid( x(:, 2), ye);

%%
A = phi(eps, sqrt( (X2-X1).^2 + (Y2-Y1).^2 ));
Ae = phi(eps, sqrt( (Xe2-Xe1).^2 + (Ye2-Ye1).^2 ));
%A = phi(eps, abs(X2 - X1));
lam = A\f;

%%

fint = Ae*lam;

err = norm(fint - fun(xe, ye), inf)/norm(fun(xe, ye), inf)

%% Plot                       
figure(2)
plot3(x(:, 1), x(:, 2), f, 'r.', 'markersize', 30), grid on, hold on
plot3(xe, ye, fint), hold off
% plot(xe, fint, 'linewidth', 3)
% plot(xe, fun(xe), '-', 'linewidth', 1), hold off
% title(['RBF interpolation with Gaussian functions and eps = ',...
%     num2str(eps, 1)], 'fontsize', 16)
% legend({'Sampled data', 'Interpolant', 'Underlying function'},...
%     'location', 'northwest', 'fontsize', 14)
figure(3)
plot3(xe, ye, fint - fun(xe, ye)), grid on

%%
% for eps = (1/2).^(0:8)
%     A = phi(eps, abs(x-x'));
%     fprintf('eps = %2.4f, condition number is %2.4e \n', eps, cond(A));
% end



