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
numIter = 10;                     % Total number of trials
grids = 30+4*(1:numIter);         % Grid sizes in one dimension
dts = 3*R^2./(2*grids.^2);        % Time step sizes           
errors = zeros(numIter, 5, 2);    % Preallocating for erros    
%% Iterations
warning('off','all')
for opOrd = 1:2
    fprintf('Iterating with Laplacian of order %1.0f.\n', 2*opOrd)
for intOrd = 1:5
    fprintf('    Iterating with interpolation order %1.0f.\n', intOrd)
    bw = rm_bandwidth(3, intOrd, opOrd/2, 1.2);     % Bandwidth
        for j = 1:numIter
            fprintf('        Iterating with grid of size %1.0f ... ', grids(j))
            goon = 1; tic
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
            try
                 Lap = nu*laplacian_3d_matrix(x, x, x, 2*opOrd, band);
            catch ME
                if (strcmp(ME.identifier,'MATLAB:sub2ind:IndexOutOfRange'))
                    errors(j, intOrd, opOrd) = NaN;
                    fprintf('Did not work \n')
                    goon = 0; tt = toc;
                end
            end
            
            if goon
                        
            % Initial Condition
            u = football(Xc, Yc, Zc); u0=u;
            % Time Integration                      
            for t = 1:Ntime
                % Explicit closest point
                unew = u + dt*Lap*u;    % Time step in embedding space 
                u = IntMat*unew;        % Extension step 
            end
            errors(j, intOrd, opOrd) = norm(u-u0*uTrueMult, inf);
            fprintf('took %4.2f second. \n', toc)
            end
        end
end
end
warning('on','all')
%

%%

errors1 = squeeze(errors(:,:,1));
errors2 = squeeze(errors(:,:,2));


figure(1)
loglog(grids, errors1, '.', 'markersize', 25)
title('Errors with Laplacian of order 2', 'fontsize', 16)
legend({'Interpolation order 1', 'Interpolation order 2',...
    'Interpolation order 3', 'Interpolation order 4',...
    'Interpolation order 5'}, 'fontsize', 14, 'location', 'southwest')


figure(2)
loglog(grids, errors2, '.', 'markersize', 25)
title('Errors with Laplacian of order 4', 'fontsize', 16)
legend({'Interpolation order 1', 'Interpolation order 2',...
    'Interpolation order 3', 'Interpolation order 4',...
    'Interpolation order 5'}, 'fontsize', 14, 'location', 'southwest')














