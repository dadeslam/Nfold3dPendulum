function vec = reorder(z)

    N = length(z)/6;
    q = z(1:length(z)/2);
    w = z(length(z)/2+1:end);
    
    vec = zeros(length(z),1);
    
    for i = 1:2*N
        if mod(i,2)==0
            vec(3*i-2:3*i) = getVec(w,i/2);
        else
            vec(3*i-2:3*i) = getVec(q,(i+1)/2);
        end
    end

end