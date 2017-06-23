%% Swift-Hohenberg Equation On a Torus    
% u_t = -Lap^2 u - Lap u + Pu + f(u), where
% f(u) = su^2 - u^3, P = -0.9, s = 1
% With stabalized Laplacian operator
% Time stepping - modified Crank-Nikolson, Adams-Bashforth
%% Computational Parameters               
N = 50;              % Number of grid points in one direction
Nt = 2e2;            % Number of time steps
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
plotgap = 2e0;       % Number of time steps between plots
bw = rm_bandwidth(3, intOrd);    % Bandwidth
%% Physical Parameters                    
R = 6;             % Major radius
r = 4;             % Minor radius
Tfinal = 1;        % Length of the simulation
P = .1;            % "Heat bath" coefficient
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
%% Parametrization                        
[xx, yy, zz] = paramTorus(128, R, r);
%% Spatial Operators                      
fprintf('Setting up operators ... '), tic
Lap     = laplacian_3d_matrix(x, x, x, opOrd, band);                    % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
%% Time Stepping Operators                
I = speye(size(Lap));               % Identity matrix 
dLap = I.*Lap;                      % Diagonal of the Laplacian 
M = (Lap - dLap)*Ext + dLap;        % Stabilized Laplace-Beltrami 
lrhs = -M^2 - M + (P-1)*I;          % Linear part of the right hand side
A = I - 9/16*dt*lrhs;               % Implicit time step operator 
[iLL, iUU] = ilu(A);                % Preconditioners for gmres       
fprintf('done after %2.2f seconds. \n', toc)
%% Initial Condition                      
u = cos(Zc).*sin(Xc); 
%u = rand(length(Zc), 1) - 0.5;  
L2norm = norm(u, 2);
%% Initial Plot                           
uplot = ExtPlot*u;
fig1 = figure(1);
lims = [-R-r R+r -R-r R+r -r r];
s = surf(xx, yy, zz, reshape(uplot, size(xx)), 'edgecolor', 'none'); 
colorbar, axis(lims), axis square, axis off, view([-37.5 70])     
title('Time = 0', 'FontSize', 20); drawnow
%% Time Integration                       
vidObj = VideoWriter('../Pictures_Movies/SH_Torus_MCNAB.avi');    % Creating video file
vidObj.FrameRate = 1;  open(vidObj);                              % Opening video file
fprintf('Staring time integration ... '), tic

% One step with Euler to get the "second IC"
[unew, ~] = gmres(I - dt*lrhs, u + dt*f(u), [], 1e-9, 100, iLL, iUU, u); 
uold = -f(u)/2 + lrhs*u/16;
u = Ext*unew;  


for t = 2:Nt       
    % Stabilized Implicit Closest Point Step
    rhs = u + dt*( 3*f(u)/2 + 3*lrhs*u/2 + uold );
    unew = gmres(A, rhs, [], 1e-9, 100, iLL, iUU, u); 
    uold = -f(u)/2 + lrhs*u/16;
    u = Ext*unew;
    
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(ExtPlot*u, size(xx));
        s.CData = uplot;  
        title(['Time = ', num2str(t*dt, 3)])
        drawnow
        writeVideo(vidObj, getframe(fig1));
    end
end
close(vidObj);  fprintf('done after %2.1f seconds \n', toc)



%