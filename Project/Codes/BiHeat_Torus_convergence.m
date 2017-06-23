%% Bi-Heat Equation On a Torus   
% with major radius R
% and minor radius r
%% Computational Parameters      
N = 46:4:80;         % Number of grid points in one direction
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
Ntr = length(N);     % Number of trials
%% Preallocation                 
errors = zeros(Ntr, 0);   % Errors
times = zeros(Ntr, 0);    % Times
%% Physical Parameters           
R = 6;            % Major radius
r = 4;            % Minor radius
L = R+r+3;        % Half of the side length of the cube
Tfinal = 1;       % Length of the simulation
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Trials                        
for j = 1:Ntr
fprintf('Starting iteration with N = %2.0f ', N(j))
% Grids                           
dx = 2*L/N(j);               % Spatial resolution
x = dx-L:dx:L;               % 1d grid
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = dx^4/36;                % Time step size
Nt = ceil(Tfinal/dt);        % Number of time steps
dt = Tfinal/Nt;         
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
fprintf('.')
% Operators                       
Lap  = laplacian_3d_matrix(x, x, x, opOrd, band);                % Laplacian
Ext  = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);        % Extension operator
dLap = speye(size(Lap)).*Lap;  M = (Lap - dLap)*Ext + dLap;      % Stabilized Laplace-Beltrami
fprintf('.')
% Initial Condition               
ICtheta = atan2(Yc, Xc);  ICphi = acos(Zc/r);  % Closest points represented in toroidal coordinates            
u = cos(3*ICphi);  utrue = u*exp(-81*Tfinal/(r^4));  L2norm = norm(u, 2); % IC and analytical solution
fprintf('.')
% Time Integration                
tic
for t = 1:Nt       
    % Ruuth-Merriman iteration
    u = u - dt*M*(M*u);     % Biharmonic
    u = Ext*u;  
    if ( norm(u, 2) > 2*L2norm )
        disp(['Instability at step ', num2str(t)])
        break
    end
end
times(j) = toc;
errors(j) = norm(u - utrue, inf);
fprintf(' done after %4.2f seconds, error is %2.4f \n', times(j), errors(j))
end
%% Saving the data               
save('../Data/Times_explicit', 'times')
save('../Data/Errors_explicit', 'errors')
%% Plots                         
figure(1)
loglog(N, errors, '.', 'markersize', 20)
title('Errors', 'fontsize', 20)

figure(2)
loglog(N, times, '.', 'markersize', 20)
title('Times', 'fontsize', 20)



%