%% Testing Structure of Interpolation Matrix        
% on circle of radius R
%% Computational parameters    
N = 30;    M = 10;   plotgap = 1;
L = 1.5;    
R = 1;
dx = 2*L/(N+1); 
intOrd = 5;
dim = 2;
%% Grids                       
x = -L:dx:L;
[X, Y] = meshgrid(x);       % Embedding grid
XX = X(:);  YY = Y(:); 
[th, rad] = cart2pol(XX, YY);
[cpx, cpy] = pol2cart(th, R);
bw = rm_bandwidth(dim, intOrd);     % Bandwidth
band = find( abs(rad - R) <= dx*bw );
cpxvec = cpx(band);   cpyvec = cpy(band);

th = linspace(0, 2*pi, 100)';
r = ones(size(th));
[xp, yp] = pol2cart(th,r);
xpv = xp(:); ypv = yp(:);
%% Operators                   
E  = interp2_matrix(x, x, cpxvec, cpyvec, intOrd, band);
%% Initial Condition and Plot
u = 5*rand(size(cpxvec));
figure(101)
plot2d_compdomain(u, X(band), Y(band), dx, dx, 101), hold on
plot(xpv, ypv, 'r-', 'linewidth', 4), hold off
title('Initial u')


%% Iterations
% for j=1:M
%     u = E*u;
%    if mod(j, plotgap) == 0
%         figure(101)
%         plot2d_compdomain(u, X(band), Y(band), dx, dx, 101), hold on
%         plot(xpv, ypv, 'r-', 'linewidth', 4), hold off
%         title(['u after ', num2str(j), ' steps'])
%         drawnow, pause
%     end
% end

%% Eigenthings
[eVects, eVals] = eig(full(E));
    
%% Iterations
vidObj = VideoWriter('Eigenvectors.avi');       % Creating video file
vidObj.FrameRate = 1;  open(vidObj);            % Opening video file
fig101 = figure(101);
for h = 1:length(band)
    plot2d_compdomain(real(eVects(:, h)), X(band), Y(band), dx, dx, 101) 
    title(['Corresponding eigenvalue: ', num2str(eVals(h, h))], 'fontsize', 16)
    drawnow, writeVideo(vidObj, getframe(fig101));
end
close(vidObj);




















    
    
    
    
    
    
    
    
    