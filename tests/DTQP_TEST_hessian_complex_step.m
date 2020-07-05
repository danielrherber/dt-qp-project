%--------------------------------------------------------------------------
% DTQP_TEST_hessian_complex_step.m
% Tests for DTQP_hessian_complex_step.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1;

% go through the tests
for k = 1:length(tests)

    clear f Df X param p t

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1 %

        % original function
        f = @(t,dummy123,in3)in3(:,1)+in3(:,3)-in3(:,4).*in3(:,2)-in3(:,5).*in3(:,2).^2.*in3(:,3);

        % exact Hessian
        D2f = cell(9,9);
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

        nt = 4000; nu = 1; ny = 2; np = 2;
        X = rand(nt,nu+3*ny+np);
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
        D2ft2 = DTQP_hessian_complex_step(f,X,t,param);
    toc

    % test analysis
    D2 = (D2ft-D2ft2); % difference
    D2(isnan(D2)) = 0;
    D2(isinf(D2)) = 0;
    disp("Max error:")
    disp(norm(D2(:),inf)) % error
    % disp(D2ft)
    % disp(' ')
    % disp(D2ft2)

end