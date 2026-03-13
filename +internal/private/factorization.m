function [q, p, r] = factorization(integer)
    factors = Inf;
    while(any(factors > 5))
        factors = factor(integer);
        integer = integer + 1;
    end
    p = sum(factors == 3);
    q = sum(factors == 2);
    r = sum(factors == 5);
end