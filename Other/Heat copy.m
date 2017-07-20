%% Heat Equation        
% on circle of radius R
% with altered extension operator
%% Computational Parameters                
N = 50;                       % Number of grid points in one direction
intOrd = 3;                   % Interpolation order
opOrd = 2;                    % Order of the spatial operator
%% Physical Parameters                     
R = 1;            % Radius
Tfinal = 5;       % Length of the simulation
nu = .1;          % Diffusivity
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Grids 
dx = 4*R/N;              % Spatial resolution
x = -2*R:dx:2*R;
[X, Y] = meshgrid(x);     % Embedding grid
dt = dx^2/4;            % Time step size
Nt = ceil(Tfinal/dt);
dt = Tfinal/Nt;         
[TH, d] = cart2pol(X, Y);                  %
band = find(abs(d - R)<=bw*dx);            % Constructing narrow band 
[Xc, Yc] = pol2cart(TH(band), R);   % Finding closest points 
%% Operators
Lap     = nu*laplacian_2d_matrix(x, x, opOrd, band);     % Laplacian
Ext     = interp2_matrix(x, x, Xc, Yc, intOrd, band);    % Extension operator
[V,D,W] = eig(full(Ext));     % Eigendecomposition
%% Initial Condition
u = cos(TH(band)); 
%% Initial Plot                
figure(1)
lims = [-1.5*R 1.5*R -1.5*R 1.5*R min(u)-0.5 max(u)+0.2];
s = scatter3(Xc(:), Yc(:), u, 50, u, 'filled'); 
colorbar, axis(lims), axis square, view([-45 30])    
title('Time = 0', 'FontSize', 20);
%% Time Integration            
tic
for t = 1:Nt              
    %u = QQ*(UU\(LL\(PP*(RR\u))));    % Time step on the embedding grid 
    u = u + dt*Ext*Lap*u;
    %u = Ext*u;   % Interpolation to the circle
    if mod(t, 10) == 0             % Plotting
        s.ZData = u;  s.CData = u;  
        title(['Time = ', num2str(t*dt)])
        drawnow
    end
end
fprintf('Time integration took %2.1f seconds \n', toc)






    

    
    
    
    
    
    
    
    
    