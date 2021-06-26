%--------------------------------------------------------------------------
% DTQP_TEST_jacobian.m
% Tests for DTQP_jacobian.m methods
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:2;
% tests = 2;

% derivative method to test
% deriv = @DTQP_jacobian_complex_step;
% deriv = @DTQP_jacobian_real_central;
deriv = @DTQP_jacobian_real_forward;

% go through the tests
for k = 1:length(tests)

    clear f Df X param auxdata t
    rng(824935)

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1 % Vanderpol (case 2, single function, no parameters)

        % original function
        f{1} = @(t,dummy123,in3)in3(:,1)+in3(:,3)-in3(:,4).*in3(:,2)-in3(:,5).*in3(:,2).^2.*in3(:,3);

        % exact Jacobian
        Df = cell(1,5);
        Df{1} = 1;
        Df{2} = @(t,dummy123,in3)-in3(:,4)-in3(:,5).*in3(:,2).*in3(:,3).*2.0;
        Df{3} = @(t,dummy123,in3)-in3(:,5).*in3(:,2).^2+1.0;
        Df{4} = @(t,dummy123,in3)-in3(:,2);
        Df{5} = @(t,dummy123,in3)-in3(:,2).^2.*in3(:,3);

        nt = 10000; nu = 1; ny = 2; np = 2;
        X = 100*rand(nt,nu+ny+np);
        param = [];
        auxdata = [];
        t = linspace(0,5,nt)';
        %------------------------------------------------------------------
        case 2 % Tumor (multiple functions, problem parameters)

        % original function
        f{1} = @(t,in2,in3)-in3(:,2).*in2(:,1).*log(in3(:,2)./in3(:,3));
        f{2} = @(t,in2,in3)-in3(:,3).*(in2(:,3)-in2(:,2)+in2(:,5).*in3(:,1)+in2(:,4).*in3(:,2).^(2.0./3.0));

        % exact Jacobian
        Df = cell(2,4);
        [Df{:}] = deal(0);
        Df{1,2} = @(t,in2,in3)-in2(:,1)-in2(:,1).*log(in3(:,2)./in3(:,3));
        Df{1,3} = @(t,in2,in3)(in3(:,2).*in2(:,1))./in3(:,3);
        Df{2,1} = @(t,in2,in3)-in2(:,5).*in3(:,3);
        Df{2,2} = @(t,in2,in3)in2(:,4).*1.0./in3(:,2).^(1.0./3.0).*in3(:,3).*(-2.0./3.0);
        Df{2,3} = @(t,in2,in3)-in2(:,3)+in2(:,2)-in2(:,5).*in3(:,1)-in2(:,4).*in3(:,2).^(2.0./3.0);

        nt = 4000; nu = 1; ny = 3; np = 0;
        X = rand(nt,nu+ny+np);

        b = 5.85; % per day
        Mew = 0.02; % per day
        G = 0.15; % per mg of dose per day
        zeta = 0.084; % per day
        D = 0.00873; % per mm^2 per day
        param = [zeta b Mew D G];
        auxdata = [];
        t = linspace(0,1.2,nt)';
        %------------------------------------------------------------------
    end

    % run the test and time
    tic
        % compute exact Jacobian
        Dfi = DTQP_QLIN_update_tmatrix(Df,[],X,param);
        Dft = DTQP_tmultiprod(Dfi,auxdata,t);
    toc

    tic
        % compute complex step Jacobian
        Dft2 = deriv(f,X,t,param);
    toc

    % test analysis
    D = Dft-Dft2; % difference
    [M,I] = max(abs(D),[],'all','linear');
    disp("Max error:")
    disp(M) % error

end