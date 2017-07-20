%% Swift-Hohenberg Equation On a Torus    
% u_t = -Lap^2 u - Lap u + Pu + f(u), where
% f(u) = su^2 - u^3, P = -0.9, s = 1
% With stabalized Laplacian operator
%% Computational Parameters               
N = 60;              % Number of grid points in one direction
Nt = 5e3;            % Number of time steps
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
plotgap = 1e2;       % Number of time steps between plots
bw = rm_bandwidth(3, intOrd);    % Bandwidth
Nplot = 100;         % Number of points in the plot in one direction
%% Physical Parameters                    
R = 6;             % Major radius
r = 4;             % Minor radius
Tfinal = 50;       % Length of the simulation
P = .9;            % "Heat bath" coefficient
s = 1;             % Coupling coefficient 
f = @(x) s*x.^2 - x.^3;          % Nonlinear part
%% Grids                                  
L = R+r+3;                   % Half of the side length of the cube
dx = 2*L/N;                  % Spatial resolution
x = dx-L:dx:L;               % 1d grid
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = Tfinal/Nt;              % Time step size
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
%% Parametrization of The Torus           
theta = linspace(0, 2*pi, Nplot);
phi = theta';
xx = (r*sin(theta) + R).*cos(phi);
yy = (r*sin(theta) + R).*sin(phi);
zz = r*cos(theta) + 0*phi;
%% Spatial Operators                      
fprintf('Setting up operators ... '), tic
Lap     = laplacian_3d_matrix(x, x, x, opOrd, band);                    % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
fI = speye(Nplot); fII = kron(fI, fI);                                  % Identity for the rectangle
Dxx = gallery('circul', sparse([-2 1 zeros(1, Nplot-3) 1]))/(dx^2);
rL = kron(fI, Dxx) + kron(Dxx, fI);                                     % Laplacian for the rectangle
%% Time Stepping Operators                
I = speye(size(Lap));               % Identity matrix 
dLap = I.*Lap;                      % Diagonal of the Laplacian 
M = (Lap - dLap)*Ext + dLap;        % Stabilized Laplace-Beltrami  
lrhs = -M^2 - M + (P-1)*I;          % Linear part of the right hand side
flrhs = -rL^2 - rL+(P-1)*fII;       % Linear part of the right hand side
At = I - dt*lrhs;                   % Implicit time step operator for the torus 
[iLLt, iUUt] = ilu(At);             % Preconditioners for gmres 
Ar = fII - dt*flrhs;                % Implicit time step operator for the torus 
[iLLr, iUUr] = ilu(Ar);             % Preconditioners for gmres 
fprintf('done after %2.2f seconds. \n', toc)                     
%% Initial Condition                      
ICtheta = atan2(Yc, Xc);        % Closest points represented
ICphi = acos(Zc/r);             % in toroidal coordinates
u = cos(2*ICphi).*sin(3*ICtheta);   % Initial condition
%% Initial Plot                           
vidObj = VideoWriter('../Pictures_Movies/SH_unwrapedTorus.avi');     % Creating video file
vidObj.FrameRate = 2;  open(vidObj);      % Opening video file
fprintf('Staring time integration ... '), tic

uplot = ExtPlot*u;   v = reshape(uplot, Nplot, Nplot);
fig1 = figure(1);
sub1 = subplot(2, 1, 1);
s1 = surf(R*phi, r*phi, reshape(uplot, Nplot, Nplot), 'edgecolor', 'none');
colorbar, axis tight, view([0 90]) 
xlabel('2\pi R', 'fontsize', 20)
ylabel('2\pi r', 'fontsize', 20)
title('Initial condition for diffusion on the torus', 'FontSize', 24);

sub2 = subplot(2, 1, 2);
s2 = surf(R*phi, r*phi, v, 'edgecolor', 'none');
colorbar, axis tight, view([0 90]) 
xlabel('2\pi R', 'fontsize', 20)
ylabel('2\pi r', 'fontsize', 20)
title('Initial condition for diffusion on the rectangle', 'FontSize', 24);
drawnow, writeVideo(vidObj, getframe(fig1));
%% Time Integration                       
v = v(:);
for t = 1:Nt       
    % Stabilized Implicit Closest Point Step
    [unew, ~] = gmres(At, u + dt*f(u), [], 1e-9, 100, iLLt, iUUt, u); 
    u = Ext*unew;  
    [v, ~] = gmres(Ar, v + dt*f(v), [], 1e-9, 100, iLLr, iUUr, v);
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(ExtPlot*u, size(xx));
        vplot = reshape(v, Nplot, Nplot);
        s1.CData = uplot;  s1.ZData = uplot;
        s2.CData = vplot;  s2.ZData = vplot;
        title(sub1, ['Solution on the torus at time = ', num2str(t*dt, 3)])
        title(sub2, ['Solution on the rectangle at time = ', num2str(t*dt, 3)])
        drawnow, writeVideo(vidObj, getframe(fig1));
    end
end
close(vidObj);  fprintf('done after %2.1f seconds \n', toc)


%