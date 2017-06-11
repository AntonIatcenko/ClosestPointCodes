
%%
n = 32;              % Space resolution
Tfinal = 1;          % Final time
Ntime = 1e3;         % Number of time steps
dt = Tfinal/Ntime;
%%

data = zeros(n, Ntime+1);

%%
x = cos(linspace(0, pi, n))';     % Chebyshev grid

u0 = (0.05./(1.05 - cos(x)));
data(:, 1) = u0;

% Compute the Dn
ep = 5;
phi = @ (ep, xi, xj) exp(-ep^2*(xi-xj).^2);
phi_xx = @ (ep, xi, xj) (-2*ep^2 + 4*ep^4*(xi-xj).^2).*exp(-ep^2*(xi-xj).^2);

[X1, X2] = meshgrid(x);

%% Operators
A = phi(ep, X2, X1);
B = phi_xx(ep, X2, X1);
D = B/A;
I = speye(n);
DI = I + dt*D;

%% BC - homogenuous Dirichlet
D(1, :) = 0;   D(end, :) = 0;

%%
for t = 1:Ntime
    
    data(:, t+1) = DI*data(:, t); 
    
end

%%
time = 0:dt:Tfinal;
figure(1)
surf(x, time, data', 'edgecolor', 'none')
colorbar


%%
%figure(4)
%plot(x, u0, '.', 'markersize', 20)




%%
%figure(5)
%Eg = eig(D);
%plot(real(Eg), imag(Eg), '.', 'markersize', 20)






