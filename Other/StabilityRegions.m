%% Stability Regions for ODE methonds       
% adaptation of Trefethens code
%
%% Adams-Bashforth                          
% s1 = subplot(2, 2, 1);
% pos1 = get(s1,'position');
% set(s1, 'Position', pos1 + [-.1 -.055 .1 .1])  
figure(101)
plot([-8 8], [0 0], 'linewidth', 2), hold on, 
plot([0 0],[-8 8], 'linewidth', 2)
z = exp(1i*pi*(0:200)/100); r = z-1;
s = 1; p1 = plot(r./s, 'linewidth', 3);                          % order 1
s = (3-1./z)/2; p2 = plot(r./s, 'linewidth', 3);                 % order 2
s = (23-16./z+5./z.^2)/12; p3 = plot(r./s, 'linewidth', 3);      % order 3
axis([-2.5 .5 -1.5 1.5]), axis square, grid on, hold off
title('Adams-Bashforth', 'fontsize', 16)
legend([p1 p2 p3], {'Order 1', 'Order 2', 'Order 3'},...
    'location', 'northwest', 'fontsize', 14)
%% Adams-Moulton                            
% s2 = subplot(2, 2, 2);
% pos2 = get(s2,'position');
% set(s2, 'Position', pos2 + [-.1 -.055 .1 .1])
figure(102)
plot([-8 8],[0 0], 'linewidth', 2), hold on 
plot([0 0],[-8 8], 'linewidth', 2)
s = (5*z+8-1./z)/12; 
p1 = plot(r./s, 'linewidth', 3);                   % order 3
s = (9*z+19-5./z+1./z.^2)/24; 
p2 = plot(r./s, 'linewidth', 3);                   % order 4 
s = (251*z+646-264./z+106./z.^2-19./z.^3)/720;
p3 = plot(r./s, 'linewidth', 3);                   % order 5
d = 1-1./z;
s = 1-d/2-d.^2/12-d.^3/24-19*d.^4/720-3*d.^5/160; 
p4 = plot(d./s, 'linewidth', 3);                   % order 6
axis([-7 1 -4 4]), axis square, grid on 
title('Adams-Moulton', 'fontsize', 16)
legend([p1 p2 p3 p4], {'Order 3', 'Order 4', 'Order 5', 'Order 6'},...
    'location', 'northwest', 'fontsize', 14)
%% Backward differentiation                 
% s3 = subplot(2, 2, 3);
% pos3 = get(s3,'position');
% set(s3, 'Position', pos3 + [-.1 -.08 .1 .1])
figure(103)
plot([-40 40],[0 0], 'linewidth', 2), hold on, 
plot([0 0],[-40 40], 'linewidth', 2) 
r = 0; 
for i = 1:6
    r = r+(d.^i)/i;     % orders 1-6
    pl(i) = plot(r, 'linewidth', 3); 
end  
axis([-15 35 -25 25]), axis square, grid on
title('backward differentiation (stable outside the curves)', 'fontsize', 16)
legend(pl, {'Order 1', 'Order 2', 'Order 3',...
    'Order 4', 'Order 5', 'Order 6'}, 'fontsize', 14)
%% Runge-Kutta                              
% s4 = subplot(2, 2, 4);
% pos4 = get(s4,'position');
% set(s4, 'Position', pos4 + [-.1 -.08 .1 .1])
figure(104)
plot([-8 8],[0 0], 'linewidth', 2), hold on
plot([0 0],[-8 8], 'linewidth', 2)
  w = 0; W = w; for i = 2:length(z)                  % order 1
    w = w-(1+w-z(i)); W = [W; w]; end
p1 = plot(W, 'linewidth', 3);
  w = 0; W = w; for i = 2:length(z)                  % order 2
    w = w-(1+w+.5*w^2-z(i)^2)/(1+w); W = [W; w];
  end
p2 = plot(W, 'linewidth', 3);
  w = 0; W = w; for i = 2:length(z)                  % order 3
    w = w-(1+w+.5*w^2+w^3/6-z(i)^3)/(1+w+w^2/2); W = [W; w];
  end
p3 = plot(W, 'linewidth', 3);
  w = 0; W = w; for i = 2:length(z)                  % order 4
    w = w-(1+w+.5*w^2+w^3/6+w.^4/24-z(i)^4)/(1+w+w^2/2+w.^3/6);
    W = [W; w]; end
p4 = plot(W, 'linewidth', 3);
  axis([-5 2 -3.5 3.5]), axis square, grid on, 
  title('Runge-Kutta', 'fontsize', 16)
  legend([p1 p2 p3 p4], {'Order 3', 'Order 4', 'Order 5', 'Order 6'},...
    'location', 'northwest', 'fontsize', 14)