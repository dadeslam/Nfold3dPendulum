function w = extractw(v)
    N = length(v)/6;
    w = zeros(3*N,1);
    for i = 1:N
        w(3*i-2:3*i) = v(6*i-2:6*i);
    end
end