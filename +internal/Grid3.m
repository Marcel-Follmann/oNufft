classdef Grid3
    properties (GetAccess = public, SetAccess = immutable)
        n
        alpha
        h

        X
        S
        gamma
        r_dyn
    end

    properties (GetAccess = private, SetAccess = immutable)
        beta
        w_vec
        q_vec
        w
    end

    methods (Access = {?oNufft.Type1, ?oNufft.Type2, ?oNufft.Type3})
        function obj = Grid3(params, x, s)
            arguments
                params oNufft.internal.Params
                x (:,1) double
                s (:,1) double
            end
            obj.X = max(abs(x(end) - x(1)));
            obj.S = max(abs(s(end) - s(1)));
            obj.n = (2*params.sigma/pi)*obj.X*obj.S + params.w;
            obj.n = 2*ceil(obj.n / 2); % round to next even integer
            obj.alpha = pi*params.w/obj.n;
            obj.h = 2*pi/obj.n;
            obj.gamma = obj.n/(2*params.sigma*obj.S);
            obj.beta = params.beta;
            obj.r_dyn = ceil(obj.calculate_r_dyn(params) +1); % Make sure that the error-bound holds!

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

        function r_dyn = calculate_r_dyn(obj, params)
            num = kernel_transform_core(0, params.beta);
            denom = kernel_transform_core(0.5*pi*params.w/params.sigma, params.beta);
            r_dyn = num / denom; 
        end
    end
end