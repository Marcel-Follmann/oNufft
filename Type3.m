classdef Type3 < handle
% TYPE3 Type-3 NUFFT operator (nonuniform -> nonuniform).
%
%   obj = oNufft.Type3(x_range, s_range)
%   obj = oNufft.Type3(..., "sigma", 2, "epsilon", 1e-6)
%
% Inputs
%   x_range - Range vector of nonuniform inputs
%   s_range - Range vector of nonuniform outputs
%
% Name-Value pairs
%   sigma   - Oversampling factor
%   epsilon - Requested accuracy
%
% Workflow
%   Scaling, Spreading and FFT are performed by apply().
%   After that, interpolate() can be used to compute the nonuniform points.
%
% Known Bugs
%   Currently the outputs of apply() don't work as intended.
%
% See also oNufft.Type3/apply, oNufft.Type3/interpolate

    properties (GetAccess = public, SetAccess = immutable)
        params oNufft.internal.Params = oNufft.internal.Params.empty
        grid oNufft.internal.Grid3 = oNufft.internal.Grid3.empty
        x_range
        s_range
        x0
        s0
    end

    properties (GetAccess = private, SetAccess = immutable)
        type_2 oNufft.Type2 = oNufft.Type2.empty
        pre_phaseing = @(x, c) c;
        post_phaseing = @(s, f) f;
    end

    properties (Access = private)
        is_ready (1,1) logical = false
    end

    methods
        function obj = Type3(x_range, s_range, options)
            arguments
                x_range (:,1) double
                s_range (:,1) double
                options.sigma (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.sigma, 2)} = 2
                options.epsilon (1,1) double = 1e-6
            end
            obj.params = oNufft.internal.Params(1, options.sigma, options.epsilon);
            obj.grid = oNufft.internal.Grid3(obj.params, x_range, s_range);
            
            obj.x_range = [x_range(1), x_range(end)];
            obj.s_range = [s_range(1), s_range(end)];
            obj.x0 = mean(obj.x_range);
            obj.s0 = mean(obj.s_range);

            if any([obj.x0, obj.s0] ~= 0)
                obj.pre_phaseing = @(x, c) obj.pre_phaseing_implementation(x, c);
                obj.post_phaseing = @(s, f) obj.post_phaseing_implementation(s, f);
            end

            obj.type_2 = oNufft.Type2(obj.grid.n, "epsilon",options.epsilon/obj.grid.r_dyn, "sigma",options.sigma);
        end

        function [s, f] = apply(obj, x, c)
            % APPLY Perform scaling, spreading and FFT stages.
            %   obj = oNufft.Type3(x_range, s_range)
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
            % See also oNufft.Type3
            arguments
                obj oNufft.Type3
                x (:,1) double
                c (:,1) double
            end
            assert(numel(x) == numel(c), "x and c must have same length");
            x = x-obj.x0;
            c = obj.pre_phaseing(x, conj(c)); % Why do i need conjugate here?

            x_scaled = (x)/obj.grid.gamma;
            [x_spreaded, c_spreaded] = obj.spread(x_scaled, c);
            [s_scaled, f_scaled] = obj.type_2.apply(x_spreaded, c_spreaded);
            f_scaled = conj(f_scaled);% Why do i need conjugate here?
            f = obj.correction(s_scaled, f_scaled);
            s = s_scaled; 

            obj.is_ready = true; 
        end

        function f = interpolate(obj, s)
            % INTERPOLATE Perform interpolation stage.
            %   obj = oNufft.Type3(x_range, s_range)
            %   obj.apply(x, c)
            %   f = obj.interpolate(s)
            %
            % Inputs
            %   s   - Nonuniform output locations
            %
            % Outputs
            %   f   - Nonuniform transform data
            %
            % See also oNufft.Type3
            arguments
                obj (1,1) oNufft.Type3
                s (:,1) double
            end
            assert(obj.is_ready, "Make sure to perform apply() before interpolation()");
            s = s - obj.s0;
            s_scaled = obj.grid.gamma*(s);
            f_scaled = obj.type_2.interpolate(obj.grid.h*s_scaled);

            f = obj.correction(s_scaled, f_scaled);
            f = obj.post_phaseing(s, f);
            % f = conj(f);
        end
    end

    methods (Access = private)
        function [x_, c_] = spread(obj, x, c)
            arguments
                obj (1,1) oNufft.Type3
                x (:,1) double
                c (:,1) double
            end
            x_ = (-floor(obj.grid.n/2):(ceil(obj.grid.n/2) - 1)).'; 
            indices = select_neigbours(x - pi, obj.params.w, obj.grid.n);
            spreaded_diracs = c .* obj.grid.periodic_kernel( x_(indices).' * obj.grid.h - x );
            indices = indices.'; % otherwise accumarray would not be working correctly
            c_ = accumarray( indices(:), spreaded_diracs(:), [obj.grid.n, 1]);
        end

        function f_ = correction(obj, s, f)
            p_k = obj.grid.h ./ obj.grid.kernel_transform(s);
            f_ = p_k .* f;
        end

        function c_ = pre_phaseing_implementation(obj, x, c)
            c_ = c.*exp(-1i*(x)*obj.s0);% Why do i need conjugate here?
        end

        function f_ = post_phaseing_implementation(obj, s, f)
            f_ = f.*exp(1i*(s + obj.s0)*obj.x0);
        end
    end
end