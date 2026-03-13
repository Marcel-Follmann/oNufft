classdef test_nufftType2 < matlab.unittest.TestCase
    properties
        type2 oNufft.Type2 = oNufft.Type2.empty
        N = 200
        epsilon = 1e-7
    end
    methods (TestClassSetup)
        function setup(testCase)
            testCase.type2 = oNufft.Type2(testCase.N, "epsilon",testCase.epsilon);
        end
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function test_frequency_random(testCase)
            k = (0:testCase.N-1).'; 
            
            M = 1000;
            x = sort(2*pi*rand(M, 1) - pi);

            n_monte_carlo = 1e3;

            f0 = NaN * ones(n_monte_carlo, 1);
            c = NaN * ones(n_monte_carlo, M);

            parfor i = 1:n_monte_carlo
                f = zeros(size(k)); 
                f0(i) = randi(floor(testCase.N/2))-1;
                f(f0(i)+1) = 1;
                
                [~] = testCase.type2.apply(k, f);
                c(i,:) = testCase.type2.interpolate(x);           
            end
            scale = 1;

            % nexttile; 
            % plot(k, f);
            % nexttile; 
            % plot(x, real(c));

            testCase.verifyEqual(real(c.'), cos(x.*f0.'), "AbsTol", scale*testCase.epsilon);
            testCase.verifyEqual(imag(c.'), sin(x.*f0.'), "AbsTol", scale*testCase.epsilon);
        end
    end

end