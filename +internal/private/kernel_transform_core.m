function phi = kernel_transform_core(z, beta)
    tmp = sqrt(beta^2 - z.^2);
    phi = (2/besseli(0, beta)) * sinh(tmp)./tmp;
end