function q = extractq(v)
    N = length(v)/6;
    q = zeros(3*N,1);
    for i = 1:N
        q(3*i-2:3*i) = v(6*i-5:6*i-3);
    end
end