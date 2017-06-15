%% Narrow band vs full embedding grid
% comparing the number of unknows with two approaches
%%
R = 1;
N = 30;                 % Number of trials
G = (1:N)*10;           % Grids
sizes = zeros(4, N);    % Preallocating

dim = 3;  % dimension
bw2 = rm_bandwidth(dim, 2);
bw3 = rm_bandwidth(dim, 3);
bw4 = rm_bandwidth(dim, 4);

%%
for j = 1:N  
    dx = 4/G(j);    x = dx-2:dx:2;
    [X, Y, Z] = meshgrid(x);
    R2 = X.^2 + Y.^2 + Z.^2; R2 = R2(:);
    sizes(1, j) = numel(X);
    sizes(2, j) = sum( R2>R-bw2 & R2<= R+bw2);
    sizes(3, j) = sum( R2>R-bw3 & R2<= R+bw3);
    sizes(4, j) = sum( R2>R-bw4 & R2<= R+bw4);
end

%%
figure(1)
plot(G, sizes', 'linewidth', 3), xlim([G(1) G(end)])
title('Number of unknows', 'fontsize', 16)
legend({'Full grid', 'Interpolation order 2', ...
    'Interpolation order 3', 'Interpolation order 4'},...
    'location', 'northwest', 'fontsize', 14) 
figure(2)
loglog(G, sizes', 'linewidth', 3), xlim([G(1) G(end)])
title('Number of unknows on log scale', 'fontsize', 16)
legend({'Full grid', 'Interpolation order 2', ...
    'Interpolation order 3', 'Interpolation order 4'},...
    'location', 'northwest', 'fontsize', 14) 

















