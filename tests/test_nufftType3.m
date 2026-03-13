classdef test_nufftType3 < matlab.unittest.TestCase
    properties
        type3 oNufft.Type3 = oNufft.Type3.empty
        AbsTol = 1e-6
        x_range = [-10, 10]
        s_range = [-pi, pi]
        epsilon = 1e-8
        n_monte_carlo = 1e3
    end
    methods (TestClassSetup)
        function setup(testCase)
            close all;
            testCase.type3 = oNufft.Type3(testCase.x_range, testCase.s_range, "epsilon",testCase.epsilon);
            testCase.AbsTol = 5*testCase.type3.grid.r_dyn * testCase.epsilon;
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function test_random_SISO(testCase)

            x = NaN * ones(testCase.n_monte_carlo,1);
            s = NaN * ones(testCase.n_monte_carlo,1);
            f = NaN * ones(testCase.n_monte_carlo,1);

            parfor i = 1:testCase.n_monte_carlo
                x(i) = diff(testCase.x_range)*rand(1)+testCase.x_range(1);
                s(i) = diff(testCase.s_range)*rand(1)+testCase.s_range(1);
                
                [~] = testCase.type3.apply(x(i), 1);
                f(i) = testCase.type3.interpolate(s(i));           
            end
            scale = 1;
            testCase.verifyEqual(real(f), cos(s.*x), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), sin(s.*x), "AbsTol", scale*testCase.AbsTol);
        end

        function test_random_SIMO(testCase)
            n_interpolation_points = 10;

            x = NaN * ones(testCase.n_monte_carlo,1);
            s = NaN * ones(testCase.n_monte_carlo,n_interpolation_points);
            f = NaN * ones(testCase.n_monte_carlo,n_interpolation_points);

            parfor i = 1:testCase.n_monte_carlo
                x(i) = diff(testCase.x_range)*rand(1)+testCase.x_range(1);
                s(i,:) = diff(testCase.s_range)*rand(n_interpolation_points,1)+testCase.s_range(1);
                
                [~] = testCase.type3.apply(x(i), 1);
                f(i,:) = testCase.type3.interpolate(s(i,:));           
            end
            scale = 1;
            testCase.verifyEqual(real(f), cos(s.*x), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), sin(s.*x), "AbsTol", scale*testCase.AbsTol);
        end

        function test_random_MISO(testCase)
            N = 10;
            x = NaN * ones(testCase.n_monte_carlo,N);
            s = NaN * ones(testCase.n_monte_carlo,1);
            f = NaN * ones(testCase.n_monte_carlo,1);

            parfor i = 1:testCase.n_monte_carlo
                x(i,:) = diff(testCase.x_range)*rand(N,1)+testCase.x_range(1);
                s(i) = diff(testCase.s_range)*rand(1)+testCase.s_range(1);
                
                [~] = testCase.type3.apply(x(i,:), ones(N,1));
                f(i) = testCase.type3.interpolate(s(i));           
            end
            signal = sum(exp(1i*x.*s),2);
            scale = N;
            testCase.verifyEqual(real(f), real(signal), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), imag(signal), "AbsTol", scale*testCase.AbsTol);
        end

        function test_random_MIMO(testCase)
            N = 5;
            n_interpolation_points = 100;

            x = NaN * ones(testCase.n_monte_carlo,N);
            s = NaN * ones(testCase.n_monte_carlo,n_interpolation_points);
            f = NaN * ones(testCase.n_monte_carlo,n_interpolation_points);

            parfor i = 1:testCase.n_monte_carlo
                x(i,:) = diff(testCase.x_range)*rand(N,1)+testCase.x_range(1);
                s(i,:) = diff(testCase.s_range)*rand(n_interpolation_points,1)+testCase.s_range(1);
                
                [~] = testCase.type3.apply(x(i,:), ones(1,N));
                f(i,:) = testCase.type3.interpolate(s(i,:));           
            end
            x = permute(x, [1, 3, 2]);
            signal = sum(exp(1i*x.*s),3);
            scale = N;

            % plot(s, real(signal), 'o'); hold on;
            % plot(s, real(f), 'x');
            testCase.verifyEqual(real(f), real(signal), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), imag(signal), "AbsTol", scale*testCase.AbsTol);
        end

        % function test_symmetric_frequencies(testCase)
        %     x = max(testCase.x_range)*[-1;1];
        %     c = 0.5*ones(2,1);
        %     s = sort(0.5*pi*rand([100,1])) + 0.5*pi;
        % 
        %     testCase.type3.apply(x, c);
        %     f = testCase.type3.interpolate(s);
        %     % plot(s, real(f)); hold on;
        %     % plot(s, imag(f));
        % end

        function test_non_centered_data(testCase)
            
            n_test_points = 100;

            x = NaN * ones(testCase.n_monte_carlo,1);
            s = NaN * ones(testCase.n_monte_carlo,n_test_points);
            f = NaN * ones(testCase.n_monte_carlo,n_test_points);
            r_dyn = NaN * ones(testCase.n_monte_carlo,1);
            parfor i = 1:testCase.n_monte_carlo
                x_range_non_centered = sort(20*rand(2,1) - 10);
                s_range_non_centered = sort(2*rand(2,1)-1);
                type3_non_centered = oNufft.Type3(x_range_non_centered, s_range_non_centered, "epsilon",testCase.epsilon);
                r_dyn(i) = type3_non_centered.grid.r_dyn;

                x(i) = mean(x_range_non_centered);
                s(i,:) = linspace(s_range_non_centered(1), s_range_non_centered(2), n_test_points).';
                [~] = type3_non_centered.apply(x(i), 1);
                f(i,:) = type3_non_centered.interpolate(s(i,:));
            end
            plot(s.', real(f).');
            scale = 1;
            testCase.verifyEqual(real(f), cos(s.*x), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), sin(s.*x), "AbsTol", scale*testCase.AbsTol);
            
        end

        function test_phase_shifted_data(testCase)
            n_test_points = 100;

            x = NaN * ones(testCase.n_monte_carlo,1);
            s = NaN * ones(testCase.n_monte_carlo,n_test_points);
            f = NaN * ones(testCase.n_monte_carlo,n_test_points);
            phase = 2*pi*rand(testCase.n_monte_carlo,1);
            c = exp(1i*phase);
            scale = norm(c, 1);
            parfor i = 1:testCase.n_monte_carlo
                x_range_non_centered = sort(20*rand(2,1) - 10);
                s_range_non_centered = sort(2*rand(2,1)-1);
                type3_non_centered = oNufft.Type3(x_range_non_centered, s_range_non_centered, "epsilon",testCase.epsilon);
    
                x(i) = mean(x_range_non_centered);
                s(i,:) = linspace(s_range_non_centered(1), s_range_non_centered(2), n_test_points).';
                
                [~] = type3_non_centered.apply(x(i), c(i));
                f(i,:) = type3_non_centered.interpolate(s(i,:));
            end
            % plot(s.', real(f).');
            testCase.verifyEqual(real(f), cos(s.*x + phase), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), sin(s.*x + phase), "AbsTol", scale*testCase.AbsTol);
        end

        function test_random_frequencies(testCase)
            n_test_points = 100;
            N = 10;

            x = NaN * ones(testCase.n_monte_carlo,N);
            s = NaN * ones(testCase.n_monte_carlo,n_test_points);
            f = NaN * ones(testCase.n_monte_carlo,n_test_points);

            parfor i = 1:testCase.n_monte_carlo
                x_range_non_centered = sort(20*rand(N,1) - 10);
                s_range_non_centered = sort(2*rand(2,1)-1);
                type3_non_centered = oNufft.Type3([x_range_non_centered(1), x_range_non_centered(end)], s_range_non_centered, "epsilon",testCase.epsilon);
    
                x(i,:) = x_range_non_centered;
                s(i,:) = linspace(s_range_non_centered(1), s_range_non_centered(2), n_test_points).';
                [~] = type3_non_centered.apply(x(i,:), ones(1, N));
                f(i,:) = type3_non_centered.interpolate(s(i,:));
            end
            x = permute(x, [1, 3, 2]);
            signal = sum(exp(1i*x.*s),3);
            scale = N;

            % nexttile; plot(s, real(signal), 'o'); hold on; plot(s, real(f), 'x');
            % nexttile; plot(s, imag(signal), 'o'); hold on; plot(s, imag(f), 'x');
            testCase.verifyEqual(real(f), real(signal), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), imag(signal), "AbsTol", scale*testCase.AbsTol);
        end

        function test_random_frequencies_with_phase_shift(testCase)
            n_test_points = 100;
            N = 10;

            x = NaN * ones(testCase.n_monte_carlo,N);
            s = NaN * ones(testCase.n_monte_carlo,n_test_points);
            f = NaN * ones(testCase.n_monte_carlo,n_test_points);
            phases = 2*pi*rand(testCase.n_monte_carlo,N);
            c = exp(1i*phases);
            scale = NaN * ones(testCase.n_monte_carlo,1);

            parfor i = 1:testCase.n_monte_carlo
                x_range_non_centered = sort(20*rand(N,1) - 10);
                s_range_non_centered = sort(2*rand(2,1)-1);
                type3_non_centered = oNufft.Type3([x_range_non_centered(1), x_range_non_centered(end)], s_range_non_centered, "epsilon",testCase.epsilon);
    
                x(i,:) = x_range_non_centered;
                s(i,:) = linspace(s_range_non_centered(1), s_range_non_centered(2), n_test_points).';
                [~] = type3_non_centered.apply(x(i,:), c(i,:));
                f(i,:) = type3_non_centered.interpolate(s(i,:));

                scale(i) = norm(c(i,:), 1);
            end
            x = permute(x, [1, 3, 2]);
            phases = permute(phases, [1, 3, 2]);
            signal = sum(exp(1i*x.*s + 1i*phases),3);
            % scale = norm(signal,1);

            % nexttile; plot(s, real(signal), 'o'); hold on; plot(s, real(f), 'x');
            % nexttile; plot(s, imag(signal), 'o'); hold on; plot(s, imag(f), 'x');
            testCase.verifyEqual(real(f), real(signal), "AbsTol", scale*testCase.AbsTol);
            testCase.verifyEqual(imag(f), imag(signal), "AbsTol", scale*testCase.AbsTol);
        end
    end

end