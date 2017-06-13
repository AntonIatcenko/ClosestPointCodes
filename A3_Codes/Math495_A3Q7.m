%% Heat Equation using Narrow Band       
%
%
%% Computational Parameters              
Nspace = 50;      % Number of grid points in one direction
Ntime = 1e3;      % Number of time steps
plotgap = 1e1;    % Number of time steps between plot 
Nplot = 128;      % Plot resolution
intOrd = 2;       % Interpolation order
opOrd = 2;        % Order of the spatial operator
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Physical Parameters                   
R = 1;            % Radius
Tfinal = 1;       % Length of the simulation
nu = 1;           % Diffusivity
%% Grids                                 
dt = Tfinal/Ntime;                % Temporal resolution
dx = 4*R/Nspace;                  % Spatial resolution
x = dx-2*R:dx:2*R;                % 1d grid
[X, Y, Z] = meshgrid(x);          % Full embedding grid
[TH, PHI, d] = cart2sph(X, Y, Z);                  %
band = find(abs(d - R)<=bw*dx);                    % Constructing narrow band 
[Xc, Yc, Zc] = sph2cart(TH(band), PHI(band), R);   % Finding closest points 
[xpl, ypl, zpl] = paramSphere(Nplot, R);           % Plotting grid
%% Operators                             
IntMat  = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);
Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);
IntPlot = interp3_matrix(x, x, x, xpl(:), ypl(:), zpl(:), intOrd, band);
%% Initial Condition and Plot            
football = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
u0 = football(Xc, Yc, Zc);  
uTrue = exp(-42*nu*Tfinal)*u0;     % Borrowed from http://bit.ly/2sdlXlM
u0plot = IntPlot*u0;

fig1 = figure(1);
s = surf(xpl, ypl, zpl, reshape(u0plot, Nplot+1, Nplot+1),...
    'edgecolor', 'none'); colorbar, axis off
title('Time 0', 'fontsize', 16)
u = u0;
%
%
%% Time Integration                      
%tic
for t = 1:Ntime

    % Explicit closest point
    unew = u + dt*Lap*u;    % Time step in embedding space 
    u = IntMat*unew;        % Extension step 

%     % RK4
%     A = Lap*u;                  Aext = IntMat*A;
%     B = Lap*(u + dt*Aext/2);    Bext = IntMat*B;
%     C = Lap*(u + dt*Bext/2);    Cext = IntMat*C;
%     D = Lap*(u + dt*Cext);      Dext = IntMat*D;
%     unew = u + dt/4*(Aext + 2*Bext + 2*Cext + D);
%     u = IntMat*unew;

    if mod(t, plotgap) == 0
        uplot = IntPlot*u;
        set(s, 'Cdata', reshape(uplot, Nplot+1, Nplot+1));
        title(gca, sprintf('Time %1.2f', dt*t), 'fontsize', 16)
        drawnow
    end

end
%toc


fprintf('Error in supremum norm - %2.1e \n', norm(u(:) - uTrue(:), inf))