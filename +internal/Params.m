classdef Params
    properties (GetAccess = public, SetAccess = immutable)
        N
        sigma
        epsilon
        w
        beta
    end

    methods (Access = {?oNufft.Type1, ?oNufft.Type2, ?oNufft.Type3})
        function obj = Params(N, sigma, epsilon)
            arguments
                N (1,1) double {mustBePositive}
                sigma (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(sigma, 2)} = 2
                epsilon (1,1) double = 1e-6
            end
            obj.N = N;
            obj.sigma = sigma;
            obj.epsilon = epsilon;
            obj.w = ceil(1.065*log10(1/epsilon))+1;
            obj.beta = 2.3 * obj.w;
        end
    end
end