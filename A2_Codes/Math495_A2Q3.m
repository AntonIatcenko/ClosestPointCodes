%% Heat equation                   
% solving heat equation on [0, Lx] x [0, Ly];
%% Computational parameters        
Nx = 20;  Ny = 20;      % Grid size
Ntime = 100;              % Number of time steps
plotgap = 1;              % Number of time steps between plots
%% Physical parameters             
Lx = 1;  Ly = 1;          % Domain size
Tfinal = 1;   nu = .1;     % Length of the simulation
dx = Lx/Nx; dy = Ly/Ny;   % Spatial resolutions
dt = Tfinal/Ntime;
x = dx:dx:Lx;  y = dy:dy:Ly;
[X, Y] = meshgrid(x, y);
%u0 = kron(x<.6, (y<.6)'); 
u0 = sin(2*pi*X).*sin(2*pi*Y);
u = u0(:);
%% Operators                       
Ix = speye(Nx);  Iy = speye(Ny); I = kron(Ix, Iy);
Dxx = (circshift(Ix, [1, 0]) - 2*Ix + circshift(Ix, [-1, 0]))/(dx^2);
Dyy = (circshift(Iy, [1, 0]) - 2*Iy + circshift(Iy, [-1, 0]))/(dy^2);        
Lap = nu*(kron(Ix, Dyy) + kron(Dxx, Iy));  [Low, Up] = lu(I - dt*Lap);
%% Initial plot                    
fig1 = figure(1);
s1 = surf(X, Y, reshape(u, Nx, Ny), 'edgecolor', 'none');
colorbar, title('Time = 0'), colormap(fig1, 'autumn')
axis([dx Lx dy Ly -1 1])%, view([0 90])
%% Time integration                
for t = 1:Ntime
    u = Up\(Low\u);
    %u = (Ix - dt/2*Dxx)\u;
    %u = (Iy - dt/2*Dyy)\(u');
    %u = u';
    
    if mod(t, plotgap) == 0
    set(s1, 'ZData', reshape(u, Nx, Ny));
    %set(s1, 'ZData', u);
    title(['Time = ', num2str(t*dt)]), drawnow
    end
end


