function vec = FuncW(z,L,m)

    q = z(1:length(z)/2);
    w = z(length(z)/2+1:end);
    
    R = assembleR(q,L,m);
    F = assembleF(q,w,m,L);
    
    vec = R\F;
end