%% Heat Equation        
% on a torus with major radius R
% and minor radius r
%% Computational Parameters                
N = 200;                       % Number of grid points in one direction
intOrd = 3;                   % Interpolation order
opOrd = 2;                    % Order of the spatial operator
%% Physical Parameters                     
R = 10;           % Major radius
r = 5;            % Minor radius
Tfinal = 1;       % Length of the simulation
nu = 1;           % Diffusivity
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Grids 
L = R+r+3;
dx = 2*L/N;                  % Spatial resolution
x = dx-L:dx:L;
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = dx^2/8;                 % Time step size
Nt = ceil(Tfinal/dt);
dt = Tfinal/Nt;         
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );
%% Parametrization
[xx, yy, zz] = paramTorus(128, R, r);
%% Operators
Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);                 % Laplacian
Ext     = interp3_matrix(x, x, x, Xc(:), Yc(:), Zc(:), intOrd, band);   % Extension operator
%Ext     = Ext(band, :);
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator
%% Initial Condition
u = cos(Zc(band)) - mean(cos(Zc(band)));  
%% Initial Plot    
uplot = ExtPlot*u;
figure(1)
lims = [-L L -L L -r-1 r+1];
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
