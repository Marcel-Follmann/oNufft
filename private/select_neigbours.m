function indices = select_neigbours(x, J, N)
    arguments
        x (1,:) double 
        J (1,1) double {mustBeInteger, mustBePositive}
        N (1,1) double {mustBeInteger, mustBePositive}
    end

    j = (1:J).';
    lb = get_left_boundary(x, J, N);
    indices = mod(lb+j, N)+1;
end

function lb = get_left_boundary(x, J, N)
    u = x*N/(2*pi);
    if mod(J,2) == 0 
        k = floor(u);
        lb = k - (J)/2; 
    else 
        k = round(u);
        lb = k - (J+1)/2; 
    end
end