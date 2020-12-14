function sol = LieEulerSE3N(vecField,action,p,h)

    k0 = zeros(length(p),1);
    k1 = vecField(k0,p);

    sol = action(exponentialSE3N(h*k1),p);   
    
end