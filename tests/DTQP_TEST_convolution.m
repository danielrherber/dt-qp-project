%--------------------------------------------------------------------------
% DTQP_TEST_convolution.m
% Test the computation of the convolution integral
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:5;
% tests = 1;

% go through the tests
for k = 1:length(tests)

    % initialize
    clear opts in A B
    in.p = [];

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        % nonsingular A, ED mesh
        opts.dt.mesh = 'ED';
        in.nt = 23;
        in.t = linspace(0,1,in.nt);
        in.h = diff(in.t);

        % constant matrices
        A = rand(3);
        B = [1 1;1 2;1 3];
        %------------------------------------------------------------------
        case 2
        % nonsingular A, USER mesh
        opts.dt.mesh = 'USER';
        in.t = [0,0.1,0.5,0.7,1];
        in.h = diff(in.t);
        in.nt = length(in.t);

        % constant matrices
        A = rand(3);
        B = [1 1;1 2;1 3];
        %------------------------------------------------------------------
        case 3
        % singular A, ED mesh
        opts.dt.mesh = 'ED';
        in.nt = 10000;
        in.t = linspace(0,1,in.nt);
        in.h = diff(in.t);

        % singular A
        A = [0 1; 0 0];
        B = [1;1];
        %------------------------------------------------------------------
        case 4
        % singular A, USER mesh
        opts.dt.mesh = 'USER';
        in.t = [0,0.1,0.5,0.7,1];
        in.h = diff(in.t);
        in.nt = length(in.t);

        % singular A
        A = [0 1; 0 0];
        B = [1;1];
        %------------------------------------------------------------------
        case 5
        % time-varying B
        opts.dt.mesh = 'LGL';
        in.nt = 100;
        in.t = DTQP_nodes_LGL(in.nt-1);
        in.h = diff(in.t);

        % time-varying matrix
        A = magic(2);
        B{1,1} = 0;
        B{2,1} = @(t) exp(-t);
        B{1,2} = @(t) sin(t);
        B{2,2} = 1;
    end

    % run the test and time
    tic
    Q = DTQP_DEFECTS_convolution(A,B,in,opts);
    toc
end