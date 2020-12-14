function [q0,w0,z0] = initializeSE3N(N)

    w0 = zeros(3*N,1);
    q0 = w0;
    z0 = rand(6*N,1);
    
    for i = 1:N
        q0(3*i-2:3*i) = rand(3,1);
        q0(3*i-2:3*i) = q0(3*i-2:3*i)/norm(q0(3*i-2:3*i),2);
        v = rand(3,1);
        w0(3*i-2:3*i) = hat(v)*q0(3*i-2:3*i);
        z0(6*i-5:6*i) = [q0(3*i-2:3*i);w0(3*i-2:3*i)];
    end

end