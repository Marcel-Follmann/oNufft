classdef Type1 < handle
% TYPE1 Type-1 NUFFT operator (nonuniform -> uniform).
%
%   obj = oNufft.Type1(N)
%   obj = oNufft.Type1(..., "sigma", 2, "epsilon", 1e-6)
%
% Inputs
%   N       - Number of output samples
%
% Name-Value pairs
%   sigma   - Oversampling factor
%   epsilon - Requested accuracy
%
% Workflow
%   The main computation is performed by apply().
%
% See also oNufft.Type1/apply

    properties (GetAccess = public, SetAccess = immutable)
        params oNufft.internal.Params = oNufft.internal.Params.empty
        grid oNufft.internal.Grid12 = oNufft.internal.Grid12.empty
    end

    methods
        function obj = Type1(N,options)
            arguments
                N (1,1) double
                options.sigma (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.sigma, 2)} = 2
                options.epsilon (1,1) double = 1e-6
            end
            obj.params = oNufft.internal.Params(N, options.sigma, options.epsilon);
            obj.grid = oNufft.internal.Grid12(obj.params); 
        end

        function [s, f] = apply(obj, x, c)
            % APPLY Perform spreading, FFT and correction stage.
            %   obj = oNufft.Type1(N)
            %   [s, f] = obj.apply(x, c)
            %
            % Inputs
            %   x       - Nonuniform sample locations
            %   c       - Sample values
            %
            % Outputs
            %   s       - Uniform output locations
            %   f       - Uniform transform data
            %
            % See also oNufft.Type1
            arguments
                obj (1,1) oNufft.Type1
                x (:,1) double
                c (:,1) double
            end
            assert(numel(x) == numel(c), "x and c must have same length");

            c_spread = obj.spread(x, c);
            [s_fft, f_fft] = obj.fft(c_spread);
            [s_trunc, f_trunc] = obj.truncate(s_fft, f_fft);
            f_corr = obj.correction(s_trunc, f_trunc);

            s = s_trunc;
            f = f_corr;
        end
    end

    methods (Access = private)
        function spread_data = spread(obj, x, c)
            arguments
                obj (1,1) oNufft.Type1
                x (:,1) double
                c (:,1) double
            end
            spread_axes = (0:obj.grid.n-1).'; 
            indices = select_neigbours(x, obj.params.w, obj.grid.n);
            spreaded_diracs = c.' .* obj.grid.periodic_kernel( spread_axes(indices) * obj.grid.h - x.' );
            spread_data = accumarray( indices(:), spreaded_diracs(:), [obj.grid.n, 1]);
        end

        function [s, f] = fft(obj, c)
            f(:,1) = fftshift( fft( c ) );
            s(:,1) = 0:obj.grid.n-1;
        end

        function [s_, f_] = truncate(obj,s, f)
            f_ = obj.grid.truncate( f );
            s_ = obj.grid.truncate( s - ( obj.grid.n / 2 ) );
        end

        function f_ = correction(obj, s, f)
            p_k = obj.grid.h ./ obj.grid.kernel_transform(s);
            f_ = p_k .* f;
        end
    end
end