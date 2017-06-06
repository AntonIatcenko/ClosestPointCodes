%% Convergence studies
M = 6;
grids = 2.^(3:M+3);
Ntime = 1; Tfinal = 1e-6;
errors = zeros(1, M);

for h=1:M
    tic, errors(h) = HeatCircle(grids(h), Ntime, Tfinal); toc
end
    
%%
loglog(errors, '.', 'markersize', 15)
