%% Heat Equation using Narrow Band         
%
%
%% Computational Parameters                
Nspace = 40;      % Number of grid points in one direction
Ntime = 1e4;      % Number of time steps
plotgap = 1e2;    % Number of time steps between plot 
Nplot = 128;      % Plot resolution
intOrd = 3;       % Interpolation order
opOrd = 2;        % Order of the spatial operator
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Physical Parameters                     
R = 1;            % Radius
Tfinal = 1;       % Length of the simulation
nu = .1;          % Diffusivity
%% Grids                                   
dt = Tfinal/Ntime;                % Temporal resolution
dx = 4*R/Nspace;                  % Spatial resolution
x = dx-2*R:dx:2*R;                % 1d grid
gam0 = 6/(dx^2);                  % Penalty parameter
[X, Y, Z] = meshgrid(x);          % Full embedding grid
[TH, PHI, d] = cart2sph(X, Y, Z);                  %
band = find(abs(d - R)<=bw*dx);                    % Constructing narrow band 
[Xc, Yc, Zc] = sph2cart(TH(band), PHI(band), R);   % Finding closest points 
[xpl, ypl, zpl] = paramSphere(Nplot, R);           % Plotting grid
%% Operators                               
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);     % Extension operator
Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);
IntPlot = interp3_matrix(x, x, x, xpl(:), ypl(:), zpl(:), intOrd, band);
Mpen = Ext*Lap - gam0*speye(length(Xc)) + gam0*Ext;
Apen = speye(length(Xc)) + dt*Mpen;
RM = speye(length(Xc)) + dt*Ext*Lap;
RMfull = Ext*RM;
%% Initial Condition                       
football = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
u0 = football(Xc, Yc, Zc);  
%% True Solution                           
uTrue = exp(-42*nu*Tfinal)*u0;     % Borrowed from http://bit.ly/2sdlXlM
u0plot = IntPlot*u0;
figure(1)
sub1 = subplot(1, 2, 1);
s1 = surf(xpl, ypl, zpl, reshape(u0plot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
title('Initial condition', 'fontsize', 16)%, colorbar
axis([-R R -R R -R R]), axis off, axis square
sub2 = subplot(1, 2, 2);
s2 = surf(xpl, ypl, zpl, reshape(u0plot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
title('Initial condition', 'fontsize', 16)%, colorbar
axis([-R R -R R -R R]), axis off, axis square
drawnow

%
%% Time Integration with Forward Euler     
u_pen = u0;  u_RM = u0;
tic
for t = 1:Ntime
    u_pen = Apen*u_pen;
    u_RM = Ext*(RM*u_RM);
    
    if mod(t, plotgap) == 0
        u_pen_plot = reshape(IntPlot*u_pen, Nplot+1, Nplot+1);
        u_RM_plot = reshape(IntPlot*u_RM, Nplot+1, Nplot+1);
        s1.CData = u_pen_plot;   s2.CData = u_RM_plot; 
        title(sub1, ['Solution with penalty method at time ', num2str(t*dt, 3)])
        title(sub2, ['Solution with Ruuth-Merriman method at time ', num2str(t*dt, 3)])
        drawnow
    end
end
t_run = toc;

%
%% Results                                 
err_pen = norm(u_pen(:) - uTrue(:), inf);
err_RM = norm(u_RM(:) - uTrue(:), inf);

fprintf('Done with time integration after %4.1f seconds, \n', t_run)
fprintf('error with penalty approach - %1.3f,\nerror with Ruuth-Merriman iteration - %1.3f \n',...
    err_pen, err_RM)






%