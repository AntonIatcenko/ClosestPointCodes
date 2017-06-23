%% Bi-Heat Equation On a Torus   
% with major radius R
% and minor radius r
%% Computational Parameters      
N = 50;             % Number of grid points in one direction
Nt = 200;            % Number of time steps
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
plotgap = 500;       % Number of time steps between plots
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
[xx, yy, zz] = paramTorus(128, R, r);
%% Operators                     
tic
Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);                 % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);            % Extension operator
ExtPlot = interp3_matrix(x, x, x, xx(:), yy(:), zz(:), intOrd, band);   % Extension operator for plotting
M = (Lap*Ext)^2;  A = speye(size(M)) + dt*M;                            % Time stepping operators
%[LL, UU, PP, QQ, RR] = lu(A);                                           % Prefactoring, takes <30 seconds
%diagA = speye(size(A)); diagA(1:length(A)+1:end) = diag(A);
[iLL, iUU] = ilu(A);     %Ltr = tril(A);
fprintf('Time to set up operators - %2.2f seconds. \n', toc)
%% Initial Condition             
u = cos(Zc) - mean(cos(Zc));    L2norm = norm(u, 2);
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
    % Implicit Closest Point Step
    %unew = A\u;  % Takes ~20 second per step
    %unew = QQ * ( UU\( LL\( PP*( RR\u ) ) ) );   % Takes ~0.5 seconds per step
    %[unew, ~,~,~,~] = gmres(A, u, 10, 1e-9, 20, [], [], u);  % Takes ~0.3 seconds per step
    unew = gmres(A, u, [], 1e-9, 100, iLL, iUU, u);  % Takes ~0.3 seconds per step
    u = Ext*unew;  
%     if ( norm(u, 2) > 2*L2norm )
%         disp(['Instability at step ', num2str(t)])
%         break
%     end
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(ExtPlot*u, size(xx));
        s.CData = uplot;  
        title(['Time = ', num2str(t*dt, 3)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)