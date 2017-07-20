%% Testing Structure of Interpolation Matrix        
% on circle of radius R
%% Computational parameters    
N = 30;   
L = 2;    
R = 1;
dx = 2*L/(N+1); 
intOrd = 3;
dim = 2;
%% Grids                       
x = -L:dx:L;
[X, Y] = meshgrid(x);       % Embedding grid
XX = X(:);  YY = Y(:); 
[th, rad] = cart2pol(X, Y);
[cpx, cpy] = pol2cart(th, R);
cpxvec = cpx(:);   cpyvec = cpy(:); 
bw = rm_bandwidth(dim, intOrd);     % Bandwidth
band = find( abs(rad - R) <= dx*bw );
cpxvec = cpxvec(band);   cpyvec = cpyvec(band);
%% Operators                   
E  = interp2_matrix(x, x, cpxvec, cpyvec, intOrd, band);
%L  = laplacian_2d_matrix(x', x', 2, 1:N^2);
%% Plots 1
figure(1)
subplot(1, 2, 1), spy(E),  title('Extension matrix', 'fontsize', 16)
subplot(1, 2, 2), spy(E'), title('Transpose of the extension matrix', 'fontsize', 16)
figure(2), imagesc(E*(E')), colorbar, title('E*E''', 'fontsize', 16)
%subplot(2, 2, 1), spy(L),  title('Laplacian matrix')
%subplot(2, 2, 2), spy(L'), title('Transpose of the Laplacian matrix')
%% SVD
[U, S, V] = svd(full(E));  % E = U*S*V', S are the singular values
%% 
figure(3), plot(1:length(E), diag(S),  '.', 'markersize', 20)
hold on, plot(1:length(E), ones(1, length(E)), 'linewidth', 1), hold off
title('Singular values', 'fontsize', 16), xlim([1 length(E)])

figure(31)
subplot(1, 2, 1), imagesc(U), colorbar, title('U')
subplot(1, 2, 2), imagesc(V), colorbar, title('V')

%%
S0 = round(S, 1);
E0 = U*S*(V');

figure(51), spy(E0)
title({'Before rounding $\tilde E$'}, 'fontsize', 16, 'interpreter', 'latex')

E1 = round(E0, 14); 
figure(52), spy(E1)
title({'Rounding $\tilde E$ to 14 digits' }, 'fontsize', 16, 'interpreter', 'latex')


figure(53), imagesc(E0 - E), colorbar
title({'Original E less unrounded $\tilde E$' }, 'fontsize', 16, 'interpreter', 'latex')


%% Eigenthings
[eVects, eVals] = eig(full(E));
rEvals = real(diag(eVals));
iEvals = imag(diag(eVals));

%%
figure(61)
plot(rEvals, iEvals, '.', 'markersize', 20)












