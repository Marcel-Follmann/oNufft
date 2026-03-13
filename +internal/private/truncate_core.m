function x = truncate_core(x, N)
    M = ceil(numel(x)/2);
    x = x(M-N/2:M+N/2-1);
end