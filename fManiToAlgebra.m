function F = fManiToAlgebra(q,w,L,m)

    R = assembleR(q,L,m);
    det(R);
    Func = assembleF(q,w,m,L);
    
    vec = @(v,i) getVec(v,i);
    
    N = length(m);
    
    V = R\Func;
    A = @(i) hat(vec(q,i))*vec(V,i);
    
    F = zeros(3*N,1);
    
    f = cell(1,N);
    
    for i = 1:N
        f{i} = A(i);
    end
    
    for i = 1:2*N
        if mod(i,2)==0
            F(3*i-2:3*i) = vec(w,i/2) ;
        else
            F(3*i-2:3*i) =f{(i+1)/2};
        end
    end
end