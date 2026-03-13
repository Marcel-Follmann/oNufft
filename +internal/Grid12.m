classdef Grid12

    properties (GetAccess = public, SetAccess = immutable)
        n
        alpha
        h
    end

    properties (GetAccess = private, SetAccess = immutable)
        N
        beta
        w_vec
        q_vec
        w
    end

    methods (Access = {?oNufft.Type1, ?oNufft.Type2, ?oNufft.Type3})
        function obj = Grid12(params)
            arguments
                params oNufft.internal.Params
            end
            obj.n = calculate_dft_size(params.sigma, params.N, params.w);
            obj.alpha = pi*params.w/obj.n;
            obj.h = 2*pi/obj.n;
            obj.beta = params.beta;
            obj.N = params.N;

            p = ceil(1.5*params.w + 2);

            [q_all, w_all] = oNufft.internal.lgwt.lgwt(2*p, -1, 1);
            obj.q_vec = q_all(1:p);
            obj.w_vec = w_all(1:p);
            obj.w = params.w;
        end

        function psi = kernel(obj, x)
            psi = kernel_core(x / obj.alpha, obj.beta);
        end


        function phi = periodic_kernel(obj, x)
            phi = periodic_kernel_core(x, obj.alpha, obj.beta);
        end

        function psi_hat = kernel_transform(obj, x)
            arguments
                obj 
                x (:,1)
            end
            psi_hat = obj.w * obj.h * sum( obj.w_vec.' .* kernel_core(obj.q_vec.', obj.beta) .* cos(obj.alpha * x .* obj.q_vec.'), 2);
        end

        function x = truncate(obj, x)
            x = truncate_core(x, obj.N);
        end
    end
end