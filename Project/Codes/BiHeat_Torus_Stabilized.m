%% Bi-Heat Equation On a Torus   
% with major radius R
% and minor radius r
% with stabalized Laplacian operator
%% Computational Parameters      
N = 50;             % Number of grid points in one direction
Nt = 100;            % Number of time steps
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
plotgap = 1;         % Number of time steps between plots
Nplot = 128;         % Number of points in the plot in one direction
%% Physical Parameters           
R = 6;            % Major radius
r = 4;            % Minor radius
Tfinal = 1;       % Length of the simulation
nu = 1;           % Diffusivity
bw = rm_bandwidth(3, intOrd);     % Bandwidth
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
[xx, yy, zz] = paramTorus(Nplot, R, r);
%% Operators                     
tic
Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);                 % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
%% Time Stepping Operators       
I = speye(size(Lap));            % Identity matrix
dLap = I.*Lap;                   % Diagonal of the Laplacian 
M = (Lap - dLap)*Ext + dLap;     % Stabilized Laplace-Beltrami    
A = I + dt*M^2;                  % Implicit time step operator   
[iLL, iUU] = ilu(A);             % Preconditioners for gmres   
fprintf('Time to set up operators - %2.2f seconds. \n', toc)
%% Initial Condition             
u = cos(Zc);    L2norm = norm(u, 2);
%% Initial Plot                  
uplot = ExtPlot*u;
figure(1)
lims = [-R-r R+r -R-r R+r -r r];
s = surf(xx, yy, zz, reshape(uplot, size(xx)), 'edgecolor', 'none'); 
colorbar, axis(lims), axis square, axis off, %view([-45 30])    
title('Time = 0', 'FontSize', 20); drawnow
%% Time Integration              
tic
for t = 1:Nt       
    % Stabilized Implicit Closest Point Step
    [unew, ~] = gmres(A, u, [], 1e-9, 100, iLL, iUU, u); 
    u = Ext*unew;  
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(ExtPlot*u, size(xx));
        s.CData = uplot;  
        title(['Time = ', num2str(t*dt, 3)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)