%% Gray-Scott Equation         
% on a rectangle [-L, L] x [-L L]
%% Computational parameters    
Nx = 100;     Ny = 100;        % Grid size
Ntime = 1e3;  plotgap = 10;    % Number of time steps
Lx = 2;  Ly = 2;               % Domain size
Tfinal = 50;                   % Length of the simulation
dx = 2*Lx/Nx; dy = 2*Ly/Ny;    % Spatial resolutions
dt = Tfinal/Ntime;             % Temporal resolution
%% Physical parameters         
eps1 = 2e-5;      eps2 = 1e-5; 
F = 0.04;         k = 0.1;
%% Grids                       
x = dx-Lx:dx:Lx;
y = dy-Ly:dy:Ly;
[X, Y] = meshgrid(x, y);      
XX = X(:);  YY = Y(:); 
%% Initial Condition           
u = 1-exp(-80*((XX+.05).^2+(YY+.02).^2)); 
v = exp(-80*((XX-.05).^2+(YY-.02).^2)); 
%% Operators                   
Ix = speye(Nx);        Iy = speye(Ny);       I = kron(Ix, Iy);
Dxx = (circshift(Ix, [1, 0]) - 2*Ix + circshift(Ix, [-1, 0]))/(dx^2);
Dyy = (circshift(Iy, [1, 0]) - 2*Iy + circshift(Iy, [-1, 0]))/(dy^2);  
Lap = kron(Ix, Dyy) + kron(Dxx, Iy); 
[Low1, Up1] = lu(I - eps1*dt*Lap);   [Low2, Up2] = lu(I - eps2*dt*Lap);
%% Initial Plot                
figure(1) 
s1 = subplot(1, 2, 1);
%[~, h1] = contourf(x, y, reshape(u, Nx, Ny), 20, 'linecolor', 'none');
h1 = surf(x, y, reshape(u, Nx, Ny), 'edgecolor', 'none'); view([0 90])
colorbar, title('u at time 0', 'Fontsize', 14), axis([dx-Lx Lx dy-Ly Ly])
s2 = subplot(1, 2, 2);
%[~, h2] = contourf(x, y, reshape(v, Nx, Ny), 20, 'linecolor', 'none');
h2 = surf(x, y, reshape(v, Nx, Ny), 'edgecolor', 'none'); view([0 90])
colorbar, title('v at time 0', 'Fontsize', 14), axis([dx-Lx Lx dy-Ly Ly])
%% Time Integration            
tic
for t = 1:Ntime         
    unew = Up1\(Low1\(u - dt*u.*(v.^2) + dt*F*(1-u)));
    vnew = Up2\(Low2\(v + dt*u.*(v.^2) - dt*(k+F)*v));
    u = unew;     v = vnew;  
    if mod(t, plotgap) == 0             % Plotting
        h1.ZData = reshape(u, Nx, Ny);
        h2.ZData = reshape(v, Nx, Ny);
        title(s1, ['u at time ', num2str(t*dt)])
        title(s2, ['v at time ', num2str(t*dt)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)