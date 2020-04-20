%--------------------------------------------------------------------------
% DTQP_TEST_tmultiprod.m
% Testing the DTQP_tmultiprod function
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

% tests = 1:8;
tests = 8;

% go through the tests
for k = 1:length(tests)

    clear A P t matrices A1 A2 A3 A4 A5

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        % single constant scalar
        A = rand(1);
        P = A;
        matrices = A;
        t = (0:4)';
        %------------------------------------------------------------------
        case 2
        % single constant matrix
        A = rand(6);
        P = A;
        matrices = A;
        t = (0:4)';
        %------------------------------------------------------------------
        case 3
        % single constant matrix (cell)
        A = rand(6);
        P = A;
        matrices = {A};
        t = (0:4)';
        %------------------------------------------------------------------
        case 4
        % single time-varying matrix
        A = @(t) 2*exp(-t);
        matrices = A;
        t = (0:4)';
        P = 2*exp(-t);
        %------------------------------------------------------------------
        case 5
        % single time-varying matrix (cell)
        A = {@(t) 2*exp(-t)};
        matrices = A;
        t = (0:4)';
        P = 2*exp(-t);
        %------------------------------------------------------------------
        case 6
        % two constant matrices
        A1 = rand(2,2);
        A2 = rand(2,2);
        P = A1*A2;
        matrices = {'prod',A1,A2};
        t = (0:4)';
        %------------------------------------------------------------------
        case 7
        % two time-varying matrices
        A1 = {@(t) sin(t),0;1,0};
        A2 = {@(t) exp(-t) + 1,1;0,0};
        matrices = {'prod',A1,A2};
        t = (0:100)';

        P = zeros(length(t),2,2);
        for i = 1:length(t)
            ti = t(i);
            P(i,:,:) = [sin(ti),0;1,0]*[exp(-ti) + 1,1;0,0];
        end
        %------------------------------------------------------------------
        case 8
        % multiple input matrices
        A1 = {@(t) sin(t),0;1,0};
        A2 = {1,1;0,0};
        A3 = rand(2,10);
        A4 = ones(10,1);
        A5 = -2;
        matrices = {'prod',A1,A2,A3,A4,A5};
        t = (0:100)';
        %------------------------------------------------------------------

    end

    % run the test and time
    tic
    A = DTQP_tmultiprod(matrices,[],t);
    toc

    % test analysis
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        assert(norm(squeeze(A(1,:,:))-P,'inf')<eps)
        %------------------------------------------------------------------
        case 2
        assert(norm(squeeze(A(1,:,:))-P,'inf')<eps)
        %------------------------------------------------------------------
        case 3
        assert(norm(squeeze(A(1,:,:))-P,'inf')<eps)
        %------------------------------------------------------------------
        case 4
        assert(norm(A-P,'inf')<eps)
        %------------------------------------------------------------------
        case 5
        assert(norm(A-P,'inf')<eps)
        %------------------------------------------------------------------
        case 6
        assert(norm(squeeze(A(1,:,:))-P,'inf')<eps)
        %------------------------------------------------------------------
        case 7
        assert(norm(A(:)-P(:),'inf')  <eps)
        %------------------------------------------------------------------
        case 8
        %------------------------------------------------------------------
    end
end