function alpha = GSS(x, d)
    
    TOL = 1e-6;
    a = (3 - sqrt(5)) / 2;
    x1 = x;
    x2 = x1 + d;
    x3 = x1 + a * (x2 - x1);
    [energy3, f3, K3] = potential(x3);
    counter = 0;
    
    while (norm(x1 - x2) > TOL)
        x4 = x3 + a * (x2 - x3);
        [energy4, f4, K4] = potential(x4);
        
        if (energy4 > energy3)
            x2 = x1;
            x1 = x4;
        else
            x1 = x3;
            x3 = x4;
        end
        energy3 = energy4;
        counter = counter + 1;
    end
    
    alpha = (x1 - x)' * d / (d' * d);
end