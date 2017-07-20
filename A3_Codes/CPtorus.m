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
%surf(x3, y3, z3, K, 'edgecolor', 'none')
surf(x3, y3, z3, K, 'edgecolor', 'none')
hold on, axis tight, axis off, view([-60 60])
title('Closest Points for a Torus', 'Fontsize', 20)

%%
a1 = -2; b1 = -7; c1 = 5;

[Pth, Pr, Pz] = cart2pol(a1, b1, c1);
[Pth1, Pr1] = cart2pol(Pr-R, Pz);
[Pr2, Pz2] = pol2cart(Pth1, r);
[cpa1, cpb1, cpc1] = pol2cart(Pth, Pr2+R, Pz2);

a2 = -5; b2 = 0; c2 = 3;

[Pth, Pr, Pz] = cart2pol(a2, b2, c2);
[Pth1, Pr1] = cart2pol(Pr-R, Pz);
[Pr2, Pz2] = pol2cart(Pth1, r);
[cpa2, cpb2, cpc2] = pol2cart(Pth, Pr2+R, Pz2);

a3 = -7; b3 = 6; c3 = 6;

[Pth, Pr, Pz] = cart2pol(a3, b3, c3);
[Pth1, Pr1] = cart2pol(Pr-R, Pz);
[Pr2, Pz2] = pol2cart(Pth1, r);
[cpa3, cpb3, cpc3] = pol2cart(Pth, Pr2+R, Pz2);

%%
plot3([a1 cpa1], [b1 cpb1], [c1 cpc1],  'r-*', 'linewidth', 2), hold on
plot3([a2 cpa2], [b2 cpb2], [c2 cpc2],  'r-*', 'linewidth', 2)
plot3([a3 cpa3], [b3 cpb3], [c3 cpc3],  'r-*', 'linewidth', 2), hold off
%plot3(cpa, cpb, cpc, 'r.', 'markersize', 20)

name = sprintf('../Project/Pictures_Movies/cpTorus.pdf');
print(name, '-dpdf', '-r0')





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
















