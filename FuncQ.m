function vec = FuncQ(z)

    q = z(1:length(z)/2);
    w = z(length(z)/2+1:end);
    
    vec = zeros(length(z)/2,1);
    
    for i = 1: length(z)/6
        vec(3*i-2:3*i) = hat(w(3*i-2:3*i))*q(3*i-2:3*i);
    end

end