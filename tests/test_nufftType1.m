classdef test_nufftType1 < matlab.unittest.TestCase
    properties
        type1 oNufft.Type1 = oNufft.Type1.empty
        epsilon = 1e-7; 
        N (1,1) double = 200;
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.type1 = oNufft.Type1(testCase.N, epsilon=testCase.epsilon);
        end
    end

    methods (TestMethodSetup)
    end

    methods (Test)
        function test_frequency_nyquist_symmetrical(testCase)
            x(:,1) = pi/2*[-1, 1];
            c(:,1) = 1/2*ones(numel(x),1);
            scale = norm(c, 1);
            [k, f_k] = testCase.type1.apply(x, c);

            % nexttile; 
            % stem(x, c);
            % nexttile; 
            % plot(k, f_k);


            testCase.verifyEqual(real(f_k), cos(k*x(end)), "AbsTol", scale*testCase.epsilon);
        end

        function test_frequency_0_Hz(testCase)
            x(:,1) = 0;
            c(:,1) = 1;
            scale = norm(c, 1);
            
            [k, f_k] = testCase.type1.apply(x, c);
            testCase.verifyEqual(real(f_k), cos(k*x(end)), "AbsTol", scale*testCase.epsilon);
        end

        function test_frequency_nyquist_random(testCase)
            n_monte_carlo = 1e3;
            c = 1/2*ones(2,1);
            scale = norm(c, 1);
            x = NaN * ones(n_monte_carlo, 2);
            k = NaN * ones(n_monte_carlo, testCase.type1.grid.n/2);
            f_k = NaN * ones(n_monte_carlo, testCase.type1.grid.n/2);
            parfor i = 1:n_monte_carlo
                x(i,:) = pi.*rand(2,1) - pi/2;
                
                
                [k(i,:), f_k(i,:)] = testCase.type1.apply(x(i,:), c);
            end

            x = permute(x, [1, 3, 2]);
            signal = sum(.5*exp(-1i*x.*k) ,3 );

            testCase.verifyEqual(real(f_k), real(signal), "AbsTol", scale*testCase.epsilon);
            testCase.verifyEqual(imag(f_k), imag(signal), "AbsTol", scale*testCase.epsilon);
        end
    end

end