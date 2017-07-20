%% Swift-Hohenberg Equation On a Rectangle
% u_t = -Lap^2 u - Lap u + Pu + f(u), where
% f(u) = su^2 - u^3, P = -0.9, s = 1
% With stabalized Laplacian operator
% Time stepping - modified SBDF2
%% Computational Parameters               
N = 50;              % Number of grid points in one direction
Nt = 1e2;            % Number of time steps
opOrd = 2;           % Order of the spatial operator
plotgap = 1e0;       % Number of time steps between plots
%% Physical Parameters                    
R = 6;             % Major radius
r = 4;             % Minor radius
Tfinal = 10;       % Length of the simulation
P = .1;            % "Heat bath" coefficient
s = 1;             % Coupling coefficient 
f = @(x) s*x.^2 - x.^3;          % Nonlinear part
%% Grids     
Lx = 2*pi*R;  Ly = 2*pi*r;   % Sides of the rectangle
L = R+r+3;                   % Half of the side length of the cube
dx = 2*L/N;                  % Spatial resolution
x = dx:dx:Lx;                % 1d grid in x direction
y = dx:dx:Ly;                % 1d grid in y direction
[X, Y] = meshgrid(x, y);     % Embedding grid
dt = Tfinal/Nt;              % Time step size
%% Spatial Operators                      
fprintf('Setting up operators ... '), tic
Ix = speye(length(x));   Iy = speye(length(y));
Dxx = gallery('circul', sparse([-2 1 zeros(1, length(x)-3) 1]))/(dx^2);
Dyy = gallery('circul', sparse([-2 1 zeros(1, length(y)-3) 1]))/(dx^2);
Lap = kron(Ix, Dyy) + kron(Dxx, Iy);
%% Time Stepping Operators                
I = speye(size(Lap));               % Identity matrix 
lrhs = -Lap^2 - Lap + (P-1)*I;      % Linear part of the right hand side
A = 3*I - 2*dt*lrhs;                % Implicit time step operator 
[iLL, iUU] = ilu(A);                % Preconditioners for gmres       
fprintf('done after %2.2f seconds. \n', toc)
%% Initial Condition                      
%u = cos(X/R).*sin(Y/r); 
u = rand(size(X)) - 0.5;  
%L2norm = norm(u, 2);
%% Initial Plot                           
fig1 = figure(1);
lims = [dx Lx dx Ly];
s = surf(X, Y, u, 'edgecolor', 'none'); 
colorbar, axis(lims), axis square%, axis off     
title('Time = 0', 'FontSize', 20); drawnow
%% Time Integration                       
vidObj = VideoWriter('../Pictures_Movies/SH_Torus_MCNAB.avi');    % Creating video file
vidObj.FrameRate = 1;  open(vidObj);                              % Opening video file
fprintf('Staring time integration ... '), tic

u = u(:);
% One step with Euler to get the "second IC"
[unew, ~] = gmres(I - dt*lrhs, u + dt*f(u), [], 1e-9, 100, iLL, iUU, u); 
uold = u; u = unew;  

% SBDF2
for t = 2:Nt       
    % Stabilized Implicit Closest Point Step
    rhs = 4*u - uold + dt*( 4*f(u) - 2*f(uold) );
    [unew, ~] = gmres(A, rhs, [], 1e-9, 100, iLL, iUU, u); 
    uold = u;  u = unew;
    
    if ( mod(t, plotgap) == 0 || ~(Nt-t) )         % Plotting
        uplot = reshape(u, size(X));
        s.CData = uplot;  s.ZData = uplot;  
        title(['Time = ', num2str(t*dt, 3)])
        drawnow
        writeVideo(vidObj, getframe(fig1));
    end
end
close(vidObj);  fprintf('done after %2.1f seconds \n', toc)



%