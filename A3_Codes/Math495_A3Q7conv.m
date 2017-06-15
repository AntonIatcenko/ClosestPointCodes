%% Heat Equation using Narrow Band       
% Convergence studies for various 
% interpolation and differentiation orders.
%% Physical Parameters                   
R = 1;            % Radius
Tfinal = 1;       % Length of the simulation
nu = 1;           % Diffusivity
%% Initial Condition                     
football = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
uTrueMult = exp(-42*nu*Tfinal);
%% Iteration Parameter                   
%intOrd = 1:5;                     % Interpolation orders
opOrd = 2;                        % Order of the spatial operator
grids = 2.^(5:6);                 % Grid sizes in one dimension
dts = 3*R^2./(2*grids.^2);        % Time step sizes
numIter = length(grids);          % Total number of trials  
errors = zeros(numIter, 5, 2);    % Preallocating for erros    
%% Iterations
warning('off','all')
for intOrd = 1:5
    disp(['Iterating with interpolation order ', num2str(intOrd, 1)])
    bw = rm_bandwidth(3, intOrd);     % Bandwidth
        for j = 1:numIter
            disp(['        Iterating with grid of size ', num2str(grids(j))])
            tic
            % Grids                                 
            dt = dts(j);                      % Temporal resolution
            Ntime = ceil(Tfinal/dt);
            dt = Tfinal/Ntime;
            dx = 3*R/grids(j);                % Spatial resolution
            x = dx-3*R/2:dx:3*R/2;            % 1d grid
            [X, Y, Z] = meshgrid(x);          % Full embedding grid
            [TH, PHI, d] = cart2sph(X, Y, Z);                  %
            band = find(abs(d - R)<=bw*dx);                    % Constructing narrow band 
            [Xc, Yc, Zc] = sph2cart(TH(band), PHI(band), R);   % Finding closest points 
            % Operators                             
             IntMat  = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);
%             try
                Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);
%             catch ME
%    statements
% end  

warning('query','last')
            
            
            % Initial Condition
            u = football(Xc, Yc, Zc); u0=u;
            % Time Integration                      
            for t = 1:Ntime
                % Explicit closest point
                unew = u + dt*Lap*u;    % Time step in embedding space 
                u = IntMat*unew;        % Extension step 
            end
            errors(j, intOrd, opOrd/2) = norm(u-u0*uTrueMult, inf);
            fprintf('        Took %4.2f second. \n', toc)
        end
end
warning('on','all')
%








