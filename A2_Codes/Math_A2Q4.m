%% Heat equation               
% on circle of radius R
%% Computational parameters    
Nx = 20;   Ny = 20;      % Grid size
Nth = 32;  dth = 2*pi/Nth;
Ntime = 1e4;               % Number of time steps
plotgap = 50;              % Number of time steps between plots
%% Physical parameters         
nu = 1;  R = 1;                % Diffusivity and Raius
Lx = 2;  Ly = 2;               % Domain size
Tfinal = 1;                    % Length of the simulation
dx = 2*Lx/Nx; dy = 2*Ly/Ny;    % Spatial resolutions
dt = Tfinal/Ntime;             % Temporal resolution
x = dx-Lx:dx:Lx;
y = dy-Ly:dy:Ly;
[X, Y] = meshgrid(x, y);       % Embedding grid
XX = X(:);  YY = Y(:);
%% Initial Condition           
IC = @(x) sin(x);  % exp(-10*x.^2);
%% Black Magic                 
theta = dth:dth:2*pi;  [Cx, Cy] = pol2cart(theta, R);
[~, ind] = min(sqrt((XX-Cx).^2 + (YY-Cy).^2), [], 2);
%% Operators                   
Ix = speye(Nx);  Iy = speye(Ny); I = kron(Ix, Iy);
Dxx = (circshift(Ix, [1, 0]) - 2*Ix + circshift(Ix, [-1, 0]))/(dx^2);
Dyy = (circshift(Iy, [1, 0]) - 2*Iy + circshift(Iy, [-1, 0]))/(dy^2);        
Lap = nu*(kron(Ix, Dyy) + kron(Dxx, Iy));
%[Low, Up] = ilu(I - dt*Lap);
%% Initial Plot                
lims = [-1.5*R 1.5*R -1.5*R 1.5*R -1.5 1];
u = IC(theta);     fig1 = figure(1);
s = scatter3(Cx, Cy, u, 50, u, 'filled'); hold on
[hh, h] = contourf(x, y, reshape(u(ind), Nx, Ny), 20, 'linecolor', 'none');
colorbar, axis square, axis(lims), axis square, view([-45 30]), hold off
h.ContourZLevel = lims(5);     title('Time = 0', 'FontSize', 20);
%% Time Integration            
for t = 1:1e3 % Ntime
    Uvec = u(ind)';                     % Closest point extension
    Uvec = (I - dt*Lap)\Uvec;           % Time step on the embedding grid
    %Uvec = Up\(Low\Uvec);
    %Uvec = Uvec + dt*Lap*Uvec;   
    Ugrid = reshape(Uvec, Nx, Ny);      % Reshaping to grid valued function 
    u = interp2(X, Y, Ugrid, Cx, Cy, 'cubic');   % Interpolation to the circle
    if mod(t, plotgap) == 0             % Plotting
        s.ZData = u;  s.CData = u;  h.ZData = Ugrid;
        title(['Time = ', num2str(t*dt)])
        drawnow
    end
end