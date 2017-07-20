function stability(methods, eigenvalues, dt)

figure(128), ax = gca;
%axis([-4 4 -4 4])
hold on

z = exp(1i*pi*(0:200)/100); r = z-1;

leg = {};

for m = 1:length(methods)
    
    cellInMethdos = string(methods(m));
    
    switch cellInMethdos
        case 'FE'       % Forward Euler
            plot(r, 'linewidth', 3)
            leg = [leg 'FE'];  
            ax.XLim = [-3 1];  ax.YLim = [-3 3];
        case 'AB1'      % Adams-Bashforth of order 1
            s = 1; plot(r./s, 'linewidth', 3)
            leg = [leg 'AB1'];
        case 'AB2'      % Adams-Bashforth of order 2
            s = (3-1./z)/2; plot(r./s, 'linewidth', 3); 
            leg = [leg 'AB2'];
        case 'AB3'      % Adams-Bashforth of order 3
            s = (23-16./z+5./z.^2)/12; plot(r./s, 'linewidth', 3); 
            leg = [leg 'AB3'];
    end
    
if nargin > 1
    if nargin == 2, dt = 1; end
    plot(dt*real(eigenvalues), dt*imag(eigenvalues), '.', 'markersize', 20)
    leg = [leg 'Eigenvalues'];
end
    
title('Stability Regions', 'fontsize', 16)
legend(leg, 'fontsize', 14)

end

plot([-4 4], [0 0], 'linewidth', 1)
plot([0 0], [-4 4], 'linewidth', 1)
hold off

end

