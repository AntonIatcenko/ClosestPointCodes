%% Gray-Scott Equation         
% on a rectangle [-L, L] x [-L L]
%% Computational parameters    
Nx = 100;     Ny = 100;        % Grid size
Ntime = 5e2;  plotgap = 5;    % Number of time steps
Lx = 2;  Ly = 2;               % Domain size
Tfinal = 5;                   % Length of the simulation
dx = 2*Lx/Nx; dy = 2*Ly/Ny;    % Spatial resolutions
dt = Tfinal/Ntime;             % Temporal resolution
%% Physical parameters         
eps1 = 2e-5;      eps2 = 1e-5; 
R = 1;    F = 0.025;     k = 0.0547;
%% Grids                       
x = dx-Lx:dx:Lx;
y = dy-Ly:dy:Ly;
[X, Y] = meshgrid(x, y);      
XX = X(:);  YY = Y(:); 
[th, ~] = cart2pol(X, Y);
[cpx, cpy] = pol2cart(th, R);
cpxvec = cpx(:);   cpyvec = cpy(:); 
%% Initial Condition           
u = exp(sin(3*th(:)));    v = cos(th(:)); 
%% Operators                   
Ix = speye(Nx);        Iy = speye(Ny);       I = kron(Ix, Iy);
Dxx = (circshift(Ix, [1, 0]) - 2*Ix + circshift(Ix, [-1, 0]))/(dx^2);
Dyy = (circshift(Iy, [1, 0]) - 2*Iy + circshift(Iy, [-1, 0]))/(dy^2);  
Lap = kron(Ix, Dyy) + kron(Dxx, Iy); 
[Low1, Up1] = lu(I - eps1*dt*Lap);   [Low2, Up2] = lu(I - eps2*dt*Lap);
%% Initial Plot                
limsU = [-1.5*R 1.5*R -1.5*R 1.5*R 0 5];
limsV = [-1.5*R 1.5*R -1.5*R 1.5*R 0 5];

fig1 = figure(5);

s1 = subplot(1, 2, 1);
circU = scatter3(cpxvec, cpyvec, u + 2, 50, u, 'filled'); hold on
%[~, shadowU] = contourf(x, y, reshape(u, Nx, Ny), 50, 'linecolor', 'none');
%shadowU = pcolor(x, y, reshape(u, Nx, Ny)); set(shadowU, 'EdgeColor', 'none');
shadowU = imagesc(x, y, reshape(u, Nx, Ny));
colorbar, axis square, axis(limsU), view([-45 30]), hold off
%shadowU.ContourZLevel = limsU(5);     
title('Time = 0', 'FontSize', 20);

s2 = subplot(1, 2, 2);
circV = scatter3(cpxvec, cpyvec, v + 2, 50, v, 'filled'); hold on
%[~, shadowV] = contourf(x, y, reshape(v, Nx, Ny), 50, 'linecolor', 'none');
%shadowV = pcolor(x, y, reshape(V, Nx, Ny)); set(shadowV, 'EdgeColor', 'none');
shadowV = imagesc(x, y, reshape(v, Nx, Ny));
colorbar, axis square, axis(limsV), view([-45 30]), hold off
%shadowV.ContourZLevel = limsV(5);     
title('Time = 0', 'FontSize', 20);
%% Time Integration            
tic
vidObj = VideoWriter(sprintf('Gray_Scott_Equation_F%1.4f_k%1.4f_.avi',...
    F, k));                               % Creating video file
vidObj.FrameRate = 10;  open(vidObj);     % Opening video file
for t = 1:Ntime         
    unew = Up1\(Low1\(u - dt*u.*(v.^2) + dt*F*(1-u)));
    vnew = Up2\(Low2\(v + dt*u.*(v.^2) - dt*(k+F)*v));
    u = interp2(X, Y, reshape(unew, Nx, Ny), cpx, cpy, 'cubic');
    v = interp2(X, Y, reshape(vnew, Nx, Ny), cpx, cpy, 'cubic');
    u = u(:);     v = v(:);  
    if mod(t, plotgap) == 0             % Plotting
        shadowU.CData = reshape(u, Nx, Ny);
        shadowV.CData = reshape(v, Nx, Ny);
        circU.ZData = u+2;  circU.CData = u;
        circV.ZData = v+2;  circV.CData = v;
        title(s1, ['u at time ', num2str(t*dt)])
        title(s2, ['v at time ', num2str(t*dt)])
        drawnow, writeVideo(vidObj, getframe(fig1));
    end
end
fprintf('Time integration took %2.1f seconds \n', toc), close(vidObj);