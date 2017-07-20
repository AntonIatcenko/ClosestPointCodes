%% Bi-Heat Equation On the Unit Sphere 
%% Computational Parameters      
N = 36:4:60;         % Number of grid points in one direction
intOrd = 4;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
Ntr = length(N);     % Number of trials
%% Preallocation                 
errors = zeros(Ntr, 0);   % Errors
times = zeros(Ntr, 0);    % Times
%% Physical Parameters           
L = 1.5;            % Half of the side length of the cube
Tfinal = .5;        % Length of the simulation
bw = rm_bandwidth(3, intOrd);       % Bandwidth
%ICfun = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
%decay = exp(-42^2*Tfinal);   


ICfun = spherefun.sphharm(2,0);     % Initial condition
decay = exp(-16*Tfinal);               % Decay of the initial data      
%% Trials                        
for j = 1:Ntr
fprintf('Starting iteration with N = %2.0f ', N(j))
% Grids                           
dx = 2*L/N(j);               % Spatial resolution
x = dx-L:dx:L;               % 1d grid
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = dx^2;                   % Time step size
Nt = ceil(Tfinal/dt);        % Number of time steps
dt = Tfinal/Nt;         
[Xc, Yc, Zc, dist] = cpSphere(X, Y, Z, 1);      % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
fprintf('.')
% Operators                       
Lap  = laplacian_3d_matrix(x, x, x, opOrd, band);                % Laplacian
Ext  = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);        % Extension operator
dLap = speye(size(Lap)).*Lap;  M = (Lap - dLap)*Ext + dLap;      % Stabilized Laplace-Beltrami
A = speye(size(Lap)) + dt*M*M;                                   % Time stepping operator
%[iUU, iLL] = ilu(A);                                             % Preconditioners for gmres
[LL, UU, PP, QQ, RR] = lu(A); 
fprintf('.')
% Initial Condition               
u = ICfun(Xc, Yc, Zc);     uTrue = decay*u;
fprintf('.')
% Time Integration                
tic
for t = 1:Nt       
    unew = QQ * ( UU\( LL\( PP*( RR\u ) ) ) );
    u = Ext*u;  
end
times(j) = toc;
errors(j) = norm(u - uTrue, inf);
fprintf(' done after %4.2f seconds, error is %2.4f \n', times(j), errors(j))
end
%% Saving the data               
save('../Data/Times_implicit', 'times')
save('../Data/Errors_implicit', 'errors')
%% Plots                         
figure(1)
loglog(N, errors, '.', 'markersize', 20)
title('Errors', 'fontsize', 20)

figure(2)
loglog(N, times, '.', 'markersize', 20)
title('Times', 'fontsize', 20)



%