%% Testing Structure of Interpolation Matrix        
% on circle of radius R
%% Computational parameters    
N = 30;   
L = 2;    
R = 1;
dx = 2*L/(N+1); 
%% Grids                       
x = -L:dx:L;
[X, Y] = meshgrid(x);       % Embedding grid
XX = X(:);  YY = Y(:); 
[th, ~] = cart2pol(X, Y);
[cpx, cpy] = pol2cart(th, R);
cpxvec = cpx(:);   cpyvec = cpy(:); 
%% Operators                   
E  = interp2_matrix(x, x, cpxvec, cpyvec);
%L  = laplacian_2d_matrix(x', x', 2, 1:N^2);
%% Plots 1
figure(1)
subplot(1, 2, 1), spy(E),  title('Extension matrix', 'fontsize', 16)
subplot(1, 2, 2), spy(E'), title('Transpose of the extension matrix', 'fontsize', 16)
%subplot(2, 2, 1), spy(L),  title('Laplacian matrix')
%subplot(2, 2, 2), spy(L'), title('Transpose of the Laplacian matrix')

%% Eigenthings
[eVects, eVals] = eig(full(E));
rEvals = real(diag(eVals));
iEvals = imag(diag(eVals));
%% SVD
[U, S, V] = svd(full(E));  % E = U*S*V', S are the singular values
%% 
figure(24), plot(1:length(E), [rEvals diag(S)],  '.', 'markersize', 20)
hold on, plot(1:length(E), ones(1, length(E))), hold off
title('.. values', 'fontsize', 16), xlim([1 length(E)])
legend({'Eigenvalues', 'Singular values'}, 'fontsize', 14)
%% Plots 2
figure(2), imagesc(E*(E')), title('E*E'''), colorbar

figure(3)
title('Eigenvalues')
plot(rEvals, iEvals, '.', 'markersize', 20)


figure(4)
subplot(1, 2, 1), spy(eVals), title('Eigenvalues')
subplot(1, 2, 2), spy(eVects), title('Eigenvectors')














