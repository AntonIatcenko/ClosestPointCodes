%% Computing Geometric Properties of a Torus
%
%
%% Parametrization
R = 3;    r = 1;
phi = @(x, y, z) x.^2 + y.^2 + z.^2 - 2*sqrt(x.^2 + y.^2) + R^2 - r^2;
%% Grids
N = 20;    dx = 10/N;
x = dx-5:dx:5;
[X, Y, Z] = meshgrid(x);

%% Operators 
I = speye(N);
D = (circshift(I, [-1, 0]) - circshift(I, [1, 0]))/(2*dx);
DX = kron(D, kron(I, I));
DY = kron(I, kron(D, I));
DZ = kron(I, kron(I, D));

%% Computing the normal
PHI = phi(X, Y, Z);

chiT = (PHI>-5).*(PHI<5);

nabPHIx = DX*(PHI(:));
nabPHIy = DY*(PHI(:));
nabPHIz = DZ*(PHI(:));



nabPHIxv = reshape(nabPHIx./nr, N, N, N);
nabPHIyv = reshape(nabPHIy./nr, N, N, N);
nabPHIzv = reshape(nabPHIz./nr, N, N, N);

nr = sqrt( nabPHIxv.^2 + nabPHIyv.^2 + nabPHIzv.^2 );


nabPHIxn = nabPHIxv./nr;
nabPHIyn = nabPHIyv./nr;
nabPHIzn = nabPHIzv./nr;


%% Parametrization 
%N = 50;
theta1 = linspace(0, 2*pi, N);
theta2 = theta1';


x3 = (r*cos(theta1) + R).*cos(theta2);
y3 = (r*cos(theta1) + R).*sin(theta2);
z3 = r*sin(theta1) + 0*theta2;


K = cos(theta1)./(R + r*cos(theta2));


%%
figure(1)
surf(x3, y3, z3, 0.*K, 'edgecolor', 'none')
colorbar, hold on


%%
figure(2)
quiver3(X, Y, Z, nabPHIxn, nabPHIyn, nabPHIzn)
axis([-5 5 -5 5 -5 5])



%%
figure(3)
quiver(X(:, :, end/2), Y(:, :, end/2), nabPHIxn(:, :, end/2), nabPHIyn(:, :, end/2))
axis([-5 5 -5 5])




















