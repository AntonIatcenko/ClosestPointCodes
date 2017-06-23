%% Spectra of Surface Biharmonic      
%
%
%% Computational Parameters           
N = 40;              % Number of grid points in one direction
intOrd = 3;          % Interpolation order
opOrd = 2;           % Order of the spatial operator
%% Physical Parameters                
R = 6;            % Major radius
r = 2;            % Minor radius
Tfinal = 2;       % Length of the simulation
nu = 1;           % Diffusivity
bw = rm_bandwidth(3, intOrd);     % Bandwidth
%% Grids                              
L = R+r+3;                   % Half of the side length of the cube
dx = 2*L/N;                  % Spatial resolution
x = dx-L:dx:L;               % 1d grid
[X, Y, Z] = meshgrid(x);     % Embedding grid
dt = dx^4/36;                % Time step size
Nt = ceil(Tfinal/dt);        % Number of time steps
dt = Tfinal/Nt;         
[Xc, Yc, Zc, dist] = cpTorus(X, Y, Z, R, r);    % Finding closest points
band = find( abs(dist) <= dx*bw );              % Constructing narrow band
Xc = Xc(band);  Yc = Yc(band);  Zc = Zc(band); 
%% Operators                          
Lap     = laplacian_3d_matrix(x, x, x, opOrd, band);            % Laplacian
Ext     = interp3_matrix(x, x, x, Xc, Yc, Zc, intOrd, band);    % Extension operator
dLap = speye(size(Lap)).*Lap;     % Diagonal of the Laplacian
M = (Lap - dLap)*Ext + dLap;      % Stabilized Laplace-Beltrami
M2 = -M^2;                        % Stabilized Biharmonic
%% Spectra of the Laplace-Beltrami
r = colamd (M);      Mfull = full(M(r, r));                  % Full M
fprintf('Size of the matrix is %4.0f by %4.0f. \n', length(M), length(M))
fprintf('Starting eigenvalue computation for Laplace-Beltrami ... '), tic
evals1 = eig(Mfull);
fprintf('took %4.4f seconds. \n', toc)
%% Plot 1
figure(1)
plot(real(evals1), imag(evals1), '.', 'markersize', 20)
title('Eigenvalues of the Laplace-Beltrami', 'fontsize', 20)
name = sprintf('../Pictures_Movies/Spectra_LapBel_Torus_%4.0f', N);
drawnow, print(name, '-dpdf', '-fillpage')
%% Spectra of the Biharmonic  
Mfull2 = full(M2);              % Full Biharmonic 
fprintf('Size of the matrix is %4.0f by %4.0f. \n', length(M), length(M))
fprintf('Starting eigenvalue computation for Biharmonic ... '), tic
evals2 = eig(Mfull2);
fprintf('took %4.4f seconds. \n', toc)
%% Plot 2                            
figure(2)
plot(real(evals2), imag(evals2), '.', 'markersize', 20)
title('Eigenvalues of the Biharmonic', 'fontsize', 20)
name = sprintf('../Pictures_Movies/Spectra_Biharm_Torus_%4.0f', N);
drawnow, print(name, '-dpdf', '-fillpage')






















