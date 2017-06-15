function stability(methods, eigenvalues)

figure(128)
plot([-8 8], [0 0], 'linewidth', 2), hold on, 
plot([0 0], [-8 8], 'linewidth', 2)

leg = {};

for m = 1:length(methods)
    
    cellInMethdos = string(methods(m));
    
    switch cellInMethdos
        case 'AB1'      % Adams-Bashforth of order 1
            z = exp(1i*pi*(0:200)/100); r = z-1;
            s = 1; plot(r./s, 'linewidth', 3)
            leg = [leg 'Adams-Bashforth of order 1'];
        case 'AB2'
            disp('AB2')
    end
title('Stability Regions')
legend(leg, 'fontsize', 14)

end
end

