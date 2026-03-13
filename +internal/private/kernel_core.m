function phi = kernel_core(z, beta)
    phi = zeros(size(z));
    idx = abs(z)<=1;

    t = 1 - z(idx).^2;
    t = max(t,0);
    phi(idx) = exp(beta*( sqrt(t) - 1 ));
end