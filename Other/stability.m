function stability(methods, eigenvalues)

figure(128)
plot([-8 8], [0 0], 'linewidth', 2), hold on, 
plot([0 0],[-8 8], 'linewidth', 2)

for m = methods

    disp(m)
    
    switch m
        case 'AB1'      % Adams-Bashforth of order 1
            z = exp(1i*pi*(0:200)/100); r = z-1;
            s = 1; plot(r./s, 'linewidth', 3)
        case 'AB2'
            disp('AB2')
    end


end
end

