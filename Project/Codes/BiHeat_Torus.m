%% Bi-Heat Equation On a Torus     
% with major radius R
% and minor radius r
%% Computational Parameters        
N = 50;              % Number of grid points in one direction
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
plotgap = 10;        % Number of time steps between plots
Nplot = 100;         % Number of points in the plot in one direction
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
dt = dx^4/36;                % Time step size
Nt = ceil(Tfinal/dt);        % Number of time steps
dt = Tfinal/Nt;         
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
%% Parametrization of The Torus    
theta = linspace(0, 2*pi, Nplot);
phi = theta';
xx = (r*sin(theta) + R).*cos(phi);
yy = (r*sin(theta) + R).*sin(phi);
zz = r*cos(theta) + 0*phi;
%% Operators                       
Lap     = laplacian_3d_matrix(x, x, x, opOrd, band);                    % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
%M = Lap*Ext;                     % Standard Laplace-Beltrami
dLap = speye(size(Lap)).*Lap;     % Diagonal of the Laplacian
M = (Lap - dLap)*Ext + dLap;      % Stabilized Laplace-Beltrami
%% Initial Condition               
ICtheta = atan2(Yc, Xc);        % Closest points represented
ICphi = acos(Zc/r);             % in toroidal coordinates
u = cos(3*ICphi);    L2norm = norm(u, 2);
ds = exp(-81*nu*plotgap*dt/(r^4));
%% Initial Plot                    
uplot = ExtPlot*u;   uflat = reshape(uplot, Nplot, Nplot);
figure(1)
sub1 = subplot(2, 1, 1);
s1 = surf(R*phi, r*phi, reshape(uplot, Nplot, Nplot), 'edgecolor', 'none');
colorbar, axis tight, view([0 90]) 
xlabel('2\pi R', 'fontsize', 20)
ylabel('2\pi r', 'fontsize', 20)
title('Initial condition for diffusion on the torus', 'FontSize', 24);

sub2 = subplot(2, 1, 2);
s2 = surf(R*phi, r*phi, uflat, 'edgecolor', 'none');
colorbar, axis tight, view([0 90]) 
xlabel('2\pi R', 'fontsize', 20)
ylabel('2\pi r', 'fontsize', 20)
title('Initial condition for diffusion on the rectangle', 'FontSize', 24);


% lims = [-R-r R+r -R-r R+r -r r];
% s = surf(xx, yy, zz, reshape(uplot, Nplot, Nplot), 'edgecolor', 'none'); 
% colorbar, axis(lims), axis square, axis off, %view([-45 30])    
%% Time Integration                
tic
for t = 1:Nt       
    % Ruuth-Merriman iteration
    u = u - dt*M*(M*u);     % Biharmonic
    %u = u + dt*M*u;         % Harmonic
    u = Ext*u;  
    if ( norm(u, 2) > 2*L2norm )
        disp(['Instability at step ', num2str(t)])
        break
    end
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        time = t*dt;
        uplot = reshape(ExtPlot*u, Nplot, Nplot);
        uflat = ds*uflat;
        s1.CData = uplot;    s1.ZData = uplot;
        title(sub1, ['Nimerical solution on the torus at time ', num2str(time, 3)])
        s2.CData = uflat;    s2.ZData = uflat;
        title(sub2, ['Analytical solution on the rectangle at time ', num2str(time, 3)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)