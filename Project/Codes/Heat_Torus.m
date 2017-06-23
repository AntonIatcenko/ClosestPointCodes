%% Heat Equation On a Torus     
% with major radius R
% and minor radius r
%% Computational Parameters     
N = 100;                      % Number of grid points in one direction
intOrd = 3;                   % Interpolation order
opOrd = 2;                    % Order of the spatial operator
%% Physical Parameters          
R = 6;            % Major radius
r = 4;            % Minor radius
Tfinal = 5;       % Length of the simulation
nu = 1;           % Diffusivity
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Grids                        
L = R+r+3;                   % Half of the side length of the cube
dx = 2*L/N;                  % Spatial resolution
x = dx-L:dx:L;               % 1d grid
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = dx^2/8;                 % Time step size
Nt = ceil(Tfinal/dt);        % Number of time steps
dt = Tfinal/Nt;         
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
%% Parametrization              
[xx, yy, zz] = paramTorus(128, R, r);
%% Operators                    
Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);                 % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
%% Initial Condition            
u = cos(Zc) - mean(cos(Zc));  
%% Initial Plot                 
uplot = ExtPlot*u;
figure(1)
lims = [-R-r R+r -R-r R+r -r r];
s = surf(xx, yy, zz, reshape(uplot, size(xx)), 'edgecolor', 'none'); 
colorbar, axis(lims), axis square, axis off, %view([-45 30])    
title('Time = 0', 'FontSize', 20);
%% Time Integration             
tic
for t = 1:Nt       
    % Ruuth-Merriman iteration
    u = u + dt*Ext*Lap*u;
    u = Ext*u;   
    if ( mod(t, 5) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(ExtPlot*u, size(xx));
        s.CData = uplot;  
        title(['Time = ', num2str(t*dt, 3)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)