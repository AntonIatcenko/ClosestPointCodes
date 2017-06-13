%% Gray-Scott Equation using Narrow Band       
%
%
%% Computational Parameters                    
Nspace = 40;      % Number of grid points in one direction
Ntime = 5e4;      % Number of time steps
plotgap = 1e3;    % Number of time steps between plot 
Nplot = 64;       % Plot resolution
intOrd = 3;       % Interpolation order
opOrd = 2;        % Order of the spatial operator
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Physical Parameters                         
R = 1.5;
F = 0.054;   k = 0.063;   nu = 1/30^2;   eta = nu/3;
f = @(u, v)  -u.*v.*v + F*(1-u);
g = @(u, v)   u.*v.*v + (F+k)*v;
%% Grids                                       
Tfinal = 50;                      % Length of the simulation
dt = Tfinal/Ntime;                % Temporal resolution
dx = 4/Nspace;                    % Spatial resolution
x = (-2:dx:2)';                   % 1d grid
[X, Y, Z] = meshgrid(x);          % Full embedding grid
[TH, PHI, d] = cart2sph(X, Y, Z);                  %
band = find(abs(d - R)<=bw*dx);                    % Constructing narrow band 
[Xc, Yc, Zc] = sph2cart(TH(band), PHI(band), R);   % Finding closest points 
[xpl, ypl, zpl] = paramSphere(Nplot, R);           % Plotting grid
%% Operators                                   
IntMat  = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);
Lap      = laplacian_3d_matrix(x', x', x', opOrd, band, band);
IntPlot = interp3_matrix(x, x, x, xpl(:), ypl(:), zpl(:), intOrd, band);
%% Initial Condition                           
football = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
%u0 = football(Xc, Yc, Zc);     v0 = football(Xc, Yc, Zc);
[u0, v0] = IC3d(Xc, Yc, Zc);
%% Initial Plot                                
u0plot = IntPlot*u0;     v0plot = IntPlot*v0;

fig1 = figure(1);   sub1 = gca;
s1 = surf(xpl, ypl, zpl, reshape(u0plot, Nplot+1, Nplot+1),...
     'edgecolor', 'none'); colorbar, axis off, axis square
% title('u at time 0', 'fontsize', 16)
% sub1 = subplot(1, 2, 1);
% s1 = surf(xpl, ypl, zpl, reshape(u0plot, Nplot+1, Nplot+1),...
%     'edgecolor', 'none'); colorbar, axis off
% title('u at time 0', 'fontsize', 16)
% sub2 = subplot(1, 2, 2);
% s2 = surf(xpl, ypl, zpl, reshape(v0plot, Nplot+1, Nplot+1),...
%     'edgecolor', 'none'); colorbar, axis off
% title('v at time 0', 'fontsize', 16)


u = u0;     v = v0;
%
%
%% Time Integration                            
%tic
for t = 1:Ntime

    % Explicit closest point
    unew = u + dt*f(u, v) + dt*nu*(Lap*u);
    vnew = v + dt*g(u, v) + dt*eta*(Lap*v);
    
    % Extension step 
    u = IntMat*unew;   v = IntMat*vnew;        

%     % RK4
%     A = Lap*u;                  Aext = IntMat*A;
%     B = Lap*(u + dt*Aext/2);    Bext = IntMat*B;
%     C = Lap*(u + dt*Bext/2);    Cext = IntMat*C;
%     D = Lap*(u + dt*Cext);      Dext = IntMat*D;
%     unew = u + dt/4*(Aext + 2*Bext + 2*Cext + D);
%     u = IntMat*unew;

    if mod(t, plotgap) == 0
        uplot = IntPlot*u;   vplot = IntPlot*v;
        set(s1, 'Cdata', reshape(uplot, Nplot+1, Nplot+1));
        %set(s2, 'Cdata', reshape(vplot, Nplot+1, Nplot+1));
        title(sub1, sprintf('u at time %1.0f', dt*t), 'fontsize', 16)
        %title(sub2, sprintf('v at time %1.2f', dt*t), 'fontsize', 16)
        drawnow
    end

end
%toc


%