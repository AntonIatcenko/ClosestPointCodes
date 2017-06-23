%% Swift-Hohenberg Equation On a Torus    
% u_t = -Lap^2 u - Lap u + Pu + f(u), where
% f(u) = su^2 - u^3, P = -0.9, s = 1
% With stabalized Laplacian operator
% Time stepping - SBDF4
%% Computational Parameters               
N = 60;              % Number of grid points in one direction
Nt = 1e2;            % Number of time steps
intOrd = 5;          % Interpolation order
opOrd = 4;           % Order of the spatial operator
plotgap = 1e0;       % Number of time steps between plots
bw = rm_bandwidth(3, intOrd);    % Bandwidth
%% Physical Parameters                    
R = 6;             % Major radius
r = 4;             % Minor radius
Tfinal = 5;        % Length of the simulation
P = .9;            % "Heat bath" coefficient
s = 1;             % Coupling coefficient 
f = @(x) s*x.^2 - x.^3;          % Nonlinear part
%% Grids                                  
L = R+r+4;                   % Half of the side length of the cube
dx = 2*L/N;                  % Spatial resolution
x = dx-L:dx:L;               % 1d grid
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = Tfinal/Nt;              % Time step size
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
%% Parametrization                        
[xx, yy, zz] = paramTorus(128, R, r);
%% Spatial Operators                      
fprintf('Setting up operators '), tic
Lap     = laplacian_3d_matrix(x, x, x, opOrd, band);                    % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
%% Time Stepping Operators                
I = speye(size(Lap));               % Identity matrix 
dLap = I.*Lap;                      % Diagonal of the Laplacian 
M = (Lap - dLap)*Ext + dLap;        % Stabilized Laplace-Beltrami 
lrhs = -M^2 - M + (P-1)*I;          % Linear part of the right hand side
A1 = I - dt^2*lrhs;                 % Implicit time step operator for Euler
A2 = 3*I - 2*dt^2*lrhs;             % Implicit time step operator for SBDF2
A4 = 25/12*I - dt*lrhs;             % Implicit time step operator for SBDF4
fprintf('.')
%[iLL1, iUU1] = ilu(A1);             % Preconditioners for gmres for Euler 
fprintf('.')
[iLL2, iUU2] = ilu(A2);             % Preconditioners for gmres for SBDF2
fprintf('.')
[iLL4, iUU4] = ilu(A4);             % Preconditioners for gmres for SBDF4 
fprintf(' done after %2.2f seconds. \n', toc)
%% Initial Condition                      
%u = cos(Zc).*sin(Xc);   
u = rand(length(Zc), 1) - 0.5;   
L2norm = norm(u, 2);
%% Initial Plot                           
uplot = ExtPlot*u;
fig1 = figure(1);
lims = [-R-r R+r -R-r R+r -r r];
s = surf(xx, yy, zz, reshape(uplot, size(xx)), 'edgecolor', 'none'); 
colorbar, axis(lims), axis square, axis off, view([-37.5 70])     
title('Time = 0', 'FontSize', 20); drawnow
%% Time Integration                       
vidObj = VideoWriter('../Pictures_Movies/SH_Torus_SBDF4.avi');    % Creating video file
vidObj.FrameRate = 1;  open(vidObj);                              % Opening video file
fprintf('Staring time integration '), tic

u_antient = u;
% Euler step to get the u(x, dt^2)
[unew, ~] = gmres(A1, u + dt^2*f(u) , [], 1e-9, 100, [], [], u); 
u = Ext*unew;  
% Stepping with SBDF2 to get the u(x, dt)
uoldtemp = u_antient;
for t = 2:floor(1/dt)
    rhs = 4*u + 4*dt^2*f(u) - uoldtemp - 2*dt^2*f(uoldtemp);
    [unew, ~] = gmres(A2, rhs, [], 1e-9, 100, iLL2, iUU2, u); 
    uoldtemp = u;    u = Ext*unew;  
end
uplot = reshape(ExtPlot*u, size(xx));
s.CData = uplot;  
title(['Time = ', num2str(dt, 3)])
drawnow, writeVideo(vidObj, getframe(fig1));
fprintf('.')
u_very_old = u;
% Stepping with SBDF2 to get the u(x, 2dt)
for t = floor(1/dt)+1:floor(2/dt)
    rhs = 4*u + 4*dt^2*f(u) - uoldtemp - 2*dt^2*f(uoldtemp);
    [unew, ~] = gmres(A2, rhs, [], 1e-9, 100, iLL2, iUU2, u); 
    uoldtemp = u;    u = Ext*unew;  
end
uplot = reshape(ExtPlot*u, size(xx));
s.CData = uplot;  
title(['Time = ', num2str(2*dt, 3)])
drawnow, writeVideo(vidObj, getframe(fig1));
fprintf('.')
u_old = u;
% Stepping with SBDF2 to get the u(x, 3dt)
for t = floor(2/dt)+1:floor(3/dt)
    rhs = 4*u + 4*dt^2*f(u) - uoldtemp - 2*dt^2*f(uoldtemp);
    [unew, ~] = gmres(A2, rhs, [], 1e-9, 100, iLL2, iUU2, u); 
    uoldtemp = u;    u = Ext*unew;  
end
uplot = reshape(ExtPlot*u, size(xx));
s.CData = uplot;  
title(['Time = ', num2str(3*dt, 3)])
drawnow, writeVideo(vidObj, getframe(fig1));
fprintf('.')
% Now we can start SBDF4
Nu_antient = dt*f(u_antient);
Nu_very_old = dt*f(u_very_old);
Nu_old = dt*f(u_old);

for t = 4:Nt       
    % Stabilized Implicit Closest Point Step
    Nu = dt*f(u);
    rhs = 4*u + 4*Nu - 3*u_old - 6*Nu_old + 4/3*u_very_old ...
        + 4*Nu_very_old -u_antient/4 - Nu_antient;
    [unew, ~] = gmres(A4, rhs, [], 1e-9, 100, iLL4, iUU4, u); 
    u_antient = u_very_old;   Nu_antient = Nu_very_old;
    u_very_old = u_old;       Nu_very_old = Nu_old;
    u_old = u;                Nu_old = Nu;
    u = Ext*unew;
    
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(ExtPlot*u, size(xx));
        s.CData = uplot;  
        title(['Time = ', num2str(t*dt, 3)])
        drawnow
        writeVideo(vidObj, getframe(fig1));
    end
end
close(vidObj);  fprintf(' done after %2.1f seconds \n', toc)



%