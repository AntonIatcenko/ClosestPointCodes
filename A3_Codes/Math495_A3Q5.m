%% Heat Equation using Narrow Band         
%
%
%% Computational Parameters                
Ns = 2.^(5:7);                 % Number of grid points in one direction
intOrd = 3;                    % Interpolation order
opOrd = 2;                     % Order of the spatial operator
numIter = length(Ns);          % Total number of trials  
errors = zeros(6, numIter);    % Preallocating for errors
times = errors;                % Preallocating for execution times
%% Physical Parameters                     
R = 1;            % Radius
Tfinal = .1;      % Length of the simulation
nu = .1;          % Diffusivity
bw = rm_bandwidth(3, intOrd);     % Bandwidth
football = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5); 
uTrueMult = exp(-42*nu*Tfinal);
%% Iterations
parfor j = 1:numIter  
    dx = 4*R/Ns(j);         % Spatial resolution
    dt = dx^2/6;            % Time step size
    Nt = ceil(Tfinal/dt);
    dt = Tfinal/Nt;         % Same for all methods
    gam = 6/(dx^2);         % Penalty parameter
    x = dx-2*R:dx:2*R;      % 1d grid
    [X, Y, Z] = meshgrid(x);          % Full embedding grid
    [TH, PHI, d] = cart2sph(X, Y, Z);                  %
    band = find(abs(d - R)<=bw*dx);                    % Constructing narrow band 
    [Xc, Yc, Zc] = sph2cart(TH(band), PHI(band), R);   % Finding closest points 
    % Operators
    Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);     % Extension operator
    Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);
    % Forward Euler             
    u0 = football(Xc, Yc, Zc);  u = u0;
    M = Ext*Lap - gam*speye(length(Xc)) + gam*Ext;
    A = speye(length(Xc)) + dt*M;
    tic
    for t = 1:Nt, u = A*u; end
    times(1, j) = toc;   
    errors(1, j) = norm(u - u0*uTrueMult, inf);
end
%% Convergence Analysis
r1e = polyfit(log(Ns), log(errors(1, :)), 1);
r2t = polyfit(log(Ns), log(times(1, :)), 1);



% %% Operators                               
% Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);     % Extension operator
% Lap     = nu*laplacian_3d_matrix(x, x, x, opOrd, band);
% IntPlot = interp3_matrix(x, x, x, xpl(:), ypl(:), zpl(:), intOrd, band);
% %% Initial Condition                       
% football = spherefun.sphharm(6,0) + sqrt(14/11)*spherefun.sphharm(6,5);
% u0 = football(Xc, Yc, Zc);  
% %% True Solution                           
% uTrue = exp(-42*nu*Tfinal)*u0;     % Borrowed from http://bit.ly/2sdlXlM
% uTruePlot = IntPlot*uTrue;
% figure(1)
% subplot(2, 3, 1)
% surf(xpl, ypl, zpl, reshape(uTruePlot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
% title('True Solution', 'fontsize', 16)%, colorbar
% axis([-R R -R R -R R]), axis off, axis square, drawnow
% 
% %
% %% Time Integration with ode45             
% disp('Computing with ode45: ')
% u = u0;
% gam = gam0;
% M = Ext*Lap - gam*speye(length(Xc)) + gam*Ext; 
% 
% tspan = [0 Tfinal];
% %opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'Stats', 'on');
% tic, [times, u_ode45] = ode45(@(u, M) M*u, tspan, u0); t45 = toc;
% e45 = norm(u(:) - uTrue(:), inf);
% 
% u_ode45Plot =  IntPlot*(u_ode45(end, :)');
% 
% figure(1)
% subplot(2, 3, 4)
% surf(xpl, ypl, zpl, reshape(u_ode45Plot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
% axis([-R R -R R -R R]), axis off, axis square%, colorbar
% title({'ode45 Solution'; sprintf('error %2.2e, time %2.2f', e45, t45)}, 'fontsize', 16)
% drawnow
% 
% %
% %% Time Integration with ode15s            
% disp('Computing with ode15s.')
% u = u0;
% gam = gam0;
% M = Ext*Lap - gam*speye(length(Xc)) + gam*Ext; 
% 
% tspan = [0 Tfinal];
% tic, [times, u_ode15s] = ode15s(@(u, M) M*u, tspan, u0); 
% t15s = toc;    e15s = norm(u(:) - uTrue(:), inf);
% 
% u_ode15sPlot =  IntPlot*(u_ode15s(end, :)');
% 
% figure(1)
% subplot(2, 3, 5)
% surf(xpl, ypl, zpl, reshape(u_ode15sPlot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
% axis([-R R -R R -R R]), axis off, axis square%, colorbar
% title({'ode15s Solution'; sprintf('error %2.2e, time %2.2f', e15s, t15s)}, 'fontsize', 16)
% drawnow
% 
% %
% %% Time Integration with ode23s            
% disp('Computing with ode23s.')
% u = u0;
% gam = gam0;
% M = Ext*Lap - gam*speye(length(Xc)) + gam*Ext; 
% opts = odeset('Jacobian', M);
% 
% tspan = [0 Tfinal];
% tic, [times, u_ode23s] = ode23s(@(u, M) M*u, tspan, u0, opts); 
% t23s = toc;    e23s = norm(u(:) - uTrue(:), inf);
% 
% u_ode23sPlot =  IntPlot*(u_ode23s(end, :)');
% 
% figure(1)
% subplot(2, 3, 6)
% surf(xpl, ypl, zpl, reshape(u_ode23sPlot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
% axis([-R R -R R -R R]), axis off, axis square%, colorbar
% title({'ode23s Solution'; sprintf('error %2.2e, time %2.2f', e23s, t23s)}, 'fontsize', 16)
% drawnow
% 
% %
% %% Time Integration "by hand"              
% disp('Computing with forward Euler.')
% u = u0;
% M = Ext*Lap - gam0*speye(length(Xc)) + gam0*Ext;
% A = speye(length(Xc)) + dt*M;
% tic
% for t = 1:Ntime, u = A*u; end
% thand = toc;   ehand = norm(u(:) - uTrue(:), inf);
% 
% u_handPlot =  IntPlot*(u);
% 
% figure(1)
% subplot(2, 3, 3)
% surf(xpl, ypl, zpl, reshape(u_handPlot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
% axis([-R R -R R -R R]), axis off, axis square%, colorbar
% title({'FE Solution'; sprintf('error %2.2e, time %2.2f', ehand, thand)}, 'fontsize', 16)
% drawnow
% 
% %
% %% Time Integration by Matrix Exponential  
% disp('Computing with matrix exponential.')
% u = u0;
% gam = gam0;
% M = Ext*Lap - gam*speye(length(Xc)) + gam*Ext; 
% 
% 
% tic, u_exp = expv(Tfinal, M, u); texp = toc;
% 
% u_expPlot =  IntPlot*u_exp;
% 
% figure(1)
% subplot(2, 3, 2)
% surf(xpl, ypl, zpl, reshape(u_ode15sPlot, Nplot+1, Nplot+1), 'edgecolor', 'none'); 
% axis([-R R -R R -R R]), axis off, axis square%, colorbar
% title({'expm Solution'; sprintf('error %2.2e, time %2.2f', e15s, t15s)}, 'fontsize', 16)
% %
% % %% Wrapper for Time-Steppers
% % function [RunTime, Err] = TimeStep(IC, Tf, sol, M, method, nt)
% % 
% %     choice = string(method);
% %     
% %     switch choice
% %         case 'FE'
% %             dt = Tf/nt;
% %             A = speye(size(M)) + dt*M;
% %             u = IC;  tic
% %             for t = 1:Ntime, u = A*u; end
% %             RunTime = toc;
% %             Err = norm(u - sol, inf);
% %         case 'BE'
% %             dt = Tf/nt;
% %             A = speye(size(M)) + dt*M;
% %             u = IC;  tic
% %             for t = 1:Ntime, u = A*u; end
% %             RunTime = toc;
% %             Err = norm(u - sol, inf);
% %     end
% % 
% % end
% 
% 
% 
% 














