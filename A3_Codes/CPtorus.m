%% Closest point for a torus
%
%% Parameters
R = 5;  r = 2;  L = R+r+1;

N = 50;
dx = 2*L/N;
x = dx-L:dx:L;

%[X, Y, Z] = meshgrid(x);



%%

theta = linspace(0, 2*pi, N);
phi = theta';


x3 = (r*cos(theta) + R).*cos(phi);
y3 = (r*cos(theta) + R).*sin(phi);
z3 = r*sin(theta) + 0*phi;

K = cos(theta)./(R + r*cos(theta'));


%%
figure(1)
surf(x3, y3, z3, K, 'edgecolor', 'none')
colorbar, hold on

%%
a = 2; b = 7; c = 5;

[Pth, Pr, Pz] = cart2pol(a, b, c);
[Pth1, Pr1] = cart2pol(Pr-R, Pz);
[Pr2, Pz2] = pol2cart(Pth1, r);
[cpa, cpb, cpc] = pol2cart(Pth, Pr2+R, Pz2);

%%
plot3(a, b, c, 'b.', 'markersize', 20)
plot3(cpa, cpb, cpc, 'r.', 'markersize', 20)


hold off




% %% Closest points
% X = x; Y = x; Z = x;
% 
% 
% [TH, RR, ZZ] = cart2pol(X, Y, Z);
% [X1, Y1, Z1] = pol2cart(TH, R, ZZ);
% [TH2, ~] = cart2pol(X1 - r , Z1);
% [Xc, Yc, Zc] = pol2cart(TH, r, Z);
% 
% %%
% [Xmesh, ~] = meshgrid(Xc, Yc);
% [Ymesh, ~] = meshgrid(Yc, Zc);
% [Zmesh, ~] = meshgrid(Zc, Xc);
















