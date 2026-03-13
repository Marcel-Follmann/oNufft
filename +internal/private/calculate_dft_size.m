function n = calculate_dft_size(sigma, N, w)
    [q, p, r] = factorization(max(sigma*N, 2*w));
    n = 2^q*3^p*5^r;
end