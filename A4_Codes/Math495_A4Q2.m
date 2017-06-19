%% Heat Equation on a Periodic Interval Using RBFs     
%
%
%% Computing the derivative of the RBF                 
syms theta thetaj eps
phi = @(eps, theta, thetaj) exp( -eps^2*(2-cos(theta - thetaj)) );
phi_xx_sym = diff(phi, theta, 2);
phi_xx = matlabFunction(phi_xx_sym);
%% Some other stuff                                    
N = 32;
Tfinal = .2;         % Final time
Ntime = 2e2;         % Number of time steps
dt = Tfinal/Ntime;
data = zeros(N, Ntime+1);
%% Computational Grid                                  
dth = 2*pi/N;
th = dth-pi:dth:pi;
%% Initial Conditions                                  
u0 = (0.05./(1.05 - cos(th)));
data(:, 1) = u0;
%% Operators                                           
ep = 5;
A = phi(ep, th, th');
B = phi_xx(ep, th, th');
D = B/A;
I = speye(N);
DI = I + dt*D/(pi^2);
%% Time Integration                                    
for t = 1:Ntime, data(:, t+1) = DI*data(:, t); end
%% Plotting                                            
time = 0:dt:Tfinal;
figure(1)
surf(th/pi, time, data', 'edgecolor', 'none')
colorbar, title('Heat equation with RBFs', 'fontsize', 16)
xlabel('Space', 'fontsize', 14)
ylabel('Time', 'fontsize', 14)
axis([th(1)/pi th(end)/pi 0 Tfinal 0 1])