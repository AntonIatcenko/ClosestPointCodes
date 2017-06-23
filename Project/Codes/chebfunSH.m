R = 9;             % Major radius
r = 6;             % Minor radius
dom = [0 2*pi*R 0 2*pi*r];
tspan = [0 10];
S = spinop2(dom, tspan);
S.lin = @(u) -2*lap(u) - biharm(u);
r = 1e-2; g = 1;
S.nonlin = @(u) (-1 + r)*u + g*u.^2 - u.^3;


% u0 = 1/20*chebfun2(@(x,y) cos(x) + sin(2*x) + sin(y) + cos(2*y), dom, 'trig');
% u0 = u0 + chebfun2(@(x,y) exp(-((x-5*pi).^2 + (y-5*pi).^2)), dom, 'trig');
% u0 = u0 + chebfun2(@(x,y) exp(-((x-5*pi).^2 + (y-15*pi).^2)), dom, 'trig');
% u0 = u0 + chebfun2(@(x,y) exp(-((x-15*pi).^2 + (y-15*pi).^2)), dom, 'trig');
% u0 = u0 + chebfun2(@(x,y) exp(-((x-15*pi).^2 + (y-5*pi).^2)), dom, 'trig');
% u0 = u0 + chebfun2(@(x,y) exp(-((x-10*pi).^2 + (y-10*pi).^2)), dom, 'trig');

u0 = chebfun2(@(x,y) cos(x/R).*sin(y/r), dom, 'trig'); 

S.init = u0;

%%
figure(1)
plot(S.init), view(0,90), colorbar %axis equal, axis off

%%
figure(2)
u = spin2(S, 96, 2e-1, 'plot', 'off');
plot(u), view(0,90), colorbar



