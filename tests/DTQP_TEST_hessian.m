%--------------------------------------------------------------------------
% DTQP_TEST_hessian.m
% Tests for DTQP_hessian.m methods
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:2;
% tests = 1;

% derivative method to test
% deriv = @DTQP_hessian_complex_step;
% deriv = @DTQP_hessian_real_central;
deriv = @DTQP_hessian_real_forward;

% go through the tests
for k = 1:length(tests)

    clear f Df X param p t
    rng(35345690)

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1 %

        % original function
        f = @(t,dummy123,in3)in3(:,1)+in3(:,3)-in3(:,4).*in3(:,2)-in3(:,5).*in3(:,2).^2.*in3(:,3);

        % exact Hessian
        D2f = cell(5,5);
        [D2f{:}] = deal(0);

        D2f{2,2} = @(t,dummy123,in3)in3(:,5).*in3(:,3).*-2.0;
        D2f{2,3} = @(t,dummy123,in3)in3(:,5).*in3(:,2).*-2.0;
        D2f{2,4} = -1;
        D2f{2,5} = @(t,dummy123,in3)in3(:,2).*in3(:,3).*-2.0;
        D2f{3,2} = @(t,dummy123,in3)in3(:,5).*in3(:,2).*-2.0;
        D2f{3,5} = @(t,dummy123,in3)-in3(:,2).^2;
        D2f{4,2} = -1;
        D2f{5,2} = @(t,dummy123,in3)in3(:,2).*in3(:,3).*-2.0;
        D2f{5,3} = @(t,dummy123,in3)-in3(:,2).^2;

        nt = 10000; nu = 1; ny = 2; np = 2;
        X = 100*rand(nt,nu+ny+np);
        param = [];
        p = [];
        t = linspace(0,5,nt)';
        %------------------------------------------------------------------
        case 2 %

        % original function
        f = @(t,dummy123,in3) sin(in3(:,1));

        % exact Hessian
        D2f = cell(1,1);
        [D2f{:}] = deal(0);

        D2f{1,1} = @(t,dummy123,in3) -sin(in3(:,1));

        nt = 4000; nu = 0; ny = 1; np = 0;
        X = rand(nt,nu+ny+np);
        param = [];
        p = [];
        t = linspace(0,5,nt)';
        %------------------------------------------------------------------
    end

    % run the test and time
    tic
        % compute exact Jacobian
        D2fi = DTQP_QLIN_update_tmatrix(D2f,[],X,param);
        D2ft = DTQP_tmultiprod(D2fi,p,t);
    toc

    tic
        % compute complex step Hessian
        D2ft2 = deriv(f,X,t,param);
    toc

    % test analysis
    D = D2ft-D2ft2; % difference
    [M,I] = max(abs(D),[],'all','linear');
    disp("Max error:")
    disp(M) % error

end