function R = assembleR(q,L,m)

    N = floor(length(q)/3);
    
    R = zeros(3*N,3*N);

    for i = 1:N
        for j = i+1:N
            R(3*i-2:3*i,3*j-2:3*j) = -sum(m(j:end)) * L(i) * L(j) * hat(q(3*i-2:3*i))*hat(q(3*j-2:3*j));
        end
    end
    
    R = R + R';
    
    for i = 1:N
        R(3*i-2:3*i , 3*i-2:3*i) = sum(m(i:end)) * L(i)^2 * eye(3);
    end
    
        

end