function Err = HeatCircle(N, Nt, Tf)               
% Solves heat equation on a circle and returns supremum of error        
nu = .1;  R = 1;  L = 2;  dx = 2*L/N;  dt = Tf/Nt;                                   
x = dx-L:dx:L;     [X, Y] = meshgrid(x, x);       
[th, ~] = cart2pol(X, Y);   [cpx, cpy] = pol2cart(th, R); 
m = 2; u = sin(m*th(:));  u_true = exp(-Tf*nu*m^2)*u; 
I = speye(N);  D = (circshift(I, [1, 0]) - 2*I + circshift(I, [-1, 0]))/(dx^2);        
Lap = nu*(kron(I, D) + kron(D, I));   [Low, Up] = lu(speye(N^2) - dt*Lap); 
for t = 1:Nt   
    u = Up\(Low\u);             % Time step on the embedding grid   
    u = interp2(X, Y, reshape(u, N, N), cpx, cpy, 'cubic');   % Interpolation to the circle
    u = u(:); 
end
Err = norm(u-u_true, inf);
end