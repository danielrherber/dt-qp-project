%--------------------------------------------------------------------------
% DTQP_TEST_jacobian_complex_step.m
% Tests for DTQP_jacobian_complex_step.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:2;
% tests = 2;

% go through the tests
for k = 1:length(tests)

    clear f Df X param p t

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1 % Vanderpol (case 2, single function, no parameters)

        % original function
        f{1} = @(t,dummy123,in3)in3(:,1)+in3(:,3)-in3(:,4).*in3(:,2)-in3(:,5).*in3(:,2).^2.*in3(:,3);

        % exact Jacobian
        Df = cell(1,9);
        Df{1} = 1;
        Df{2} = @(t,dummy123,in3)-in3(:,4)-in3(:,5).*in3(:,2).*in3(:,3).*2.0;
        Df{3} = @(t,dummy123,in3)-in3(:,5).*in3(:,2).^2+1.0;
        Df{4} = @(t,dummy123,in3)-in3(:,2);
        Df{5} = @(t,dummy123,in3)-in3(:,2).^2.*in3(:,3);
        [Df{6:9}] = deal(0);

        nt = 10; nu = 1; ny = 2; np = 2;
        X = rand(nt,nu+3*ny+np);
        param = [];
        p = [];
        t = linspace(0,5,nt)';
        %------------------------------------------------------------------
        case 2 % Tumor (multiple functions, problem parameters)

        % original function
        f{1} = @(t,in2,in3)-in3(:,2).*in2(:,1).*log(in3(:,2)./in3(:,3));
        f{2} = @(t,in2,in3)-in3(:,3).*(in2(:,3)-in2(:,2)+in2(:,5).*in3(:,1)+in2(:,4).*in3(:,2).^(2.0./3.0));

        % exact Jacobian
        Df = cell(2,10);
        [Df{:}] = deal(0);
        Df{1,2} = @(t,in2,in3)-in2(:,1)-in2(:,1).*log(in3(:,2)./in3(:,3));
        Df{1,3} = @(t,in2,in3)(in3(:,2).*in2(:,1))./in3(:,3);
        Df{2,1} = @(t,in2,in3)-in2(:,5).*in3(:,3);
        Df{2,2} = @(t,in2,in3)in2(:,4).*1.0./in3(:,2).^(1.0./3.0).*in3(:,3).*(-2.0./3.0);
        Df{2,3} = @(t,in2,in3)-in2(:,3)+in2(:,2)-in2(:,5).*in3(:,1)-in2(:,4).*in3(:,2).^(2.0./3.0);

        nt = 4000; nu = 1; ny = 3; np = 0;
        X = rand(nt,nu+3*ny+np);

        b = 5.85; % per day
        Mew = 0.02; % per day
        G = 0.15; % per mg of dose per day
        zeta = 0.084; % per day
        D = 0.00873; % per mm^2 per day
        param = [zeta b Mew D G];
        p = [];
        t = linspace(0,1.2,nt)';
        %------------------------------------------------------------------
    end

    % run the test and time
    tic
        % compute exact Jacobian
        Dfi = DTQP_QLIN_update_tmatrix(Df,[],X,param);
        Dft = DTQP_tmultiprod(Dfi,p,t);
    toc

    tic
        % compute complex step Jacobian
        Dft2 = DTQP_jacobian_complex_step(f,X,t,param);
    toc

    % test analysis
    D = Dft-Dft2; % difference
    disp("Max error:")
    disp(norm(D(:),inf)) % error
    % disp(Dft)
    % disp(Dft2)

end