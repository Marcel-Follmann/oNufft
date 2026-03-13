function phi = periodic_kernel_core(z, alpha, beta)
    x = z - 2*pi*round(z/(2*pi));
    phi = sum(kernel_core(x / alpha, beta),3);
end