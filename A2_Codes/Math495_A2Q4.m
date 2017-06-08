%% Heat equation               
% on circle of radius R
%% Computational parameters    
Nx = 100;     Ny = 100;        % Grid size
Ntime = 1e2;  plotgap = 1;    % Number of time steps
%% Physical parameters         
nu = 1;  R = 1;                % Diffusivity and Raius
Lx = 2;  Ly = 2;               % Domain size
Tfinal = 1;                    % Length of the simulation
dx = 2*Lx/Nx; dy = 2*Ly/Ny;    % Spatial resolutions
dt = Tfinal/Ntime;             % Temporal resolution
%% Grids                       
x = dx-Lx:dx:Lx;
y = dy-Ly:dy:Ly;
[X, Y] = meshgrid(x, y);       % Embedding grid
XX = X(:);  YY = Y(:); 
[th, ~] = cart2pol(X, Y);
[cpx, cpy] = pol2cart(th, R);
cpxvec = cpx(:);   cpyvec = cpy(:); 
%% Initial Condition           
IC = @(x) exp(-2*x.^2);  % sin(4*x);
%% Operators                   
Ix = speye(Nx);  Iy = speye(Ny); I = kron(Ix, Iy);
Dxx = (circshift(Ix, [1, 0]) - 2*Ix + circshift(Ix, [-1, 0]))/(dx^2);
Dyy = (circshift(Iy, [1, 0]) - 2*Iy + circshift(Iy, [-1, 0]))/(dy^2);  
Lap = nu*(kron(Ix, Dyy) + kron(Dxx, Iy)); [Low, Up] = lu(I - dt*Lap);
%% Initial Plot                
u = IC(th(:));     figure(1)
lims = [-1.5*R 1.5*R -1.5*R 1.5*R min(u)-0.5 max(u)+0.2];
s = scatter3(cpxvec, cpyvec, u, 50, u, 'filled'); hold on
[~, h] = contourf(x, y, IC(atan2(Y, X)), 20, 'linecolor', 'none');
colorbar, axis square, axis(lims), axis square, view([-45 30]), hold off
h.ContourZLevel = lims(5);     title('Time = 0', 'FontSize', 20);
%% Time Integration            
tic
for t = 1:Ntime         
    u = Up\(Low\u);             % Time step on the embedding grid   
    u = interp2(X, Y, reshape(u, Nx, Ny), cpx, cpy, 'cubic');   % Interpolation to the circle
    u = u(:);
    if mod(t, plotgap) == 0             % Plotting
        s.ZData = u;  s.CData = u;  
        h.ZData = reshape(u, Nx, Ny);
        title(['Time = ', num2str(t*dt)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)