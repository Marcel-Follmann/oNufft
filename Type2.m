classdef Type2 < handle
% TYPE2 Type-2 NUFFT operator (uniform -> nonuniform).
%
%   obj = oNufft.Type2(N)
%   obj = oNufft.Type2(..., "sigma", 2, "epsilon", 1e-6)
%
% Inputs
%   N       - Number of output samples
%
% Name-Value pairs
%   sigma   - Oversampling factor
%   epsilon - Requested accuracy
%
% Workflow
%   Scaling and FFT are performed by apply().
%   After that, interpolate() can be used to compute the nonuniform points.
%
% Known Bugs
%   Currently only works for even N.
%
% See also oNufft.Type2/apply, oNufft.Type2/interpolate
    properties (GetAccess = public, SetAccess = immutable)
        params oNufft.internal.Params = oNufft.internal.Params.empty
        grid oNufft.internal.Grid12 = oNufft.internal.Grid12.empty
    end

    properties (Access = private)
        f_nyquist
        s_nyquist

        is_ready (1,1) logical = false
    end

    methods
        function obj = Type2(N, options)
            arguments
                N (:,1) double
                options.sigma (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.sigma, 2)} = 2
                options.epsilon (1,1) double = 1e-6
            end
            obj.params = oNufft.internal.Params(N, options.sigma, options.epsilon);
            obj.grid = oNufft.internal.Grid12(obj.params);
        end

        function [s, f] = apply(obj, x, c)
            % APPLY Perform scaling and FFT stages.
            %   obj = oNufft.Type2(N)
            %   [s, f] = obj.apply(x, c)
            %
            % Inputs
            %   x       - Uniform sample locations
            %   c       - Sample values
            %
            % Outputs
            %   s       - Uniform output locations
            %   f       - Uniform transform data
            %
            % See also oNufft.Type2
            arguments
                obj oNufft.Type2
                x (:,1) double
                c (:,1) double
            end
            assert(numel(x) == numel(c), "x and c must have same length");

            % Prescale
            c_scaled = obj.prescale(x, c);
            
            if x(1) < 0 
                c_scaled = obj.symetric_padding(c_scaled);
            end

            [obj.s_nyquist, obj.f_nyquist] = obj.fft(c_scaled);
            
            
            s = obj.s_nyquist*obj.grid.h;
            f = obj.f_nyquist;
            obj.is_ready = true;
        end
        
        function f = interpolate(obj, s)
            % INTERPOLATE Perform interpolation stage.
            %   obj = oNufft.Type2(N)
            %   obj.apply(x, c)
            %   f = obj.interpolate(s)
            %
            % Inputs
            %   s   - Nonuniform output locations
            %
            % Outputs
            %   f   - Nonuniform transform data
            %
            % See also oNufft.Type2
            arguments
                obj (1,1) oNufft.Type2 
                s (:,1) double
            end
            if not(obj.is_ready)
                error("oNufft.Type2 not ready: Need to 'apply' before interpolation is possible");
            end

            indices = select_neigbours(s - pi, ...
                obj.params.w, ...
                obj.grid.n);

            window = obj.grid.periodic_kernel( ...
                obj.s_nyquist(indices)*obj.grid.h - s.' ...
                );

            f(:,1) = sum( ...
                obj.f_nyquist(indices) .* window , ...
                1);
            
            f = conj(f); % Add selection for forward / backward fft.
            
        end
    end

    methods(Access = private)
        function c_scaled = prescale(obj, x, c)
            c_scaled = obj.grid.h * c./ obj.grid.kernel_transform(x);
        end

        function c_padded = symetric_padding(obj, c)
            d = obj.grid.n-obj.params.N;
            padL = floor(d/2);
            c_padded = ...
                [zeros(padL, 1); c; zeros(d - padL, 1)];
            c_padded = ifftshift(c_padded);
        end

        function [s, f] = fft(obj, c)
            f(:,1) = fftshift( ...
                fft(c, obj.grid.n) ...
                );

            s = ...
                ( -floor(obj.grid.n/2):( ceil(obj.grid.n/2) - 1 ) ).';
        end
    end
end