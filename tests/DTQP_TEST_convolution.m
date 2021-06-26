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

tests = 1:9;
% tests = 1;

% go through the tests
for k = 1:length(tests)

    % initialize
    clear opts in A B
    in.auxdata = [];

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
        %------------------------------------------------------------------
        case 6
        % time-varying B with ny > nu
        opts.dt.mesh = 'ED';
        in.nt = 100;
        in.t = linspace(0,1,in.nt);
        in.h = diff(in.t);

        % time-varying matrix
        A = magic(3);
        B{1,1} = 0;
        B{2,1} = @(t) exp(-t);
        B{3,1} = @(t) sin(t);
        B{1,2} = 1;
        B{2,2} = 0;
        B{3,2} = @(t) exp(-t);
        %------------------------------------------------------------------
        case 7
        % time-varying B with ny < nu
        opts.dt.mesh = 'ED';
        in.nt = 10;
        in.t = linspace(0,1,in.nt);
        in.h = diff(in.t);

        % time-varying matrix
        A = magic(2);
        B{1,1} = 0;
        B{2,1} = @(t) exp(-t);
        B{1,2} = 1;
        B{2,2} = 0;
        B{1,3} = @(t) sin(t);
        B{2,3} = @(t) exp(-t);
        %------------------------------------------------------------------
        case 8
        % constant B with ny < nu
        opts.dt.mesh = 'ED';
        in.nt = 10;
        in.t = linspace(0,1,in.nt);
        in.h = diff(in.t);

        A = magic(2);
        % cell form
        B{1,1} = 0;
        B{2,1} = 1;
        B{1,2} = 1;
        B{2,2} = 0;
        B{1,3} = 1;
        B{2,3} = 1;

        % double form
        % B = [0 1 1;1 0 1]; % should be the same as the cell form above
        %------------------------------------------------------------------
        case 9
        % constant B with ny < nu
        opts.dt.mesh = 'ED';
        in.nt = 10;
        in.t = linspace(0,1,in.nt);
        in.h = diff(in.t);

        A = magic(2);
        % product of matrices
        B1 = eye(2);
        B2{1,1} = 0;
        B2{2,1} = 1;
        B2{1,2} = 1;
        B2{2,2} = 0;
        B2{1,3} = @(t) sin(t);
        B2{2,3} = @(t) exp(-t);
        B{1} = 'prod';
        B{2} = B1;
        B{3} = B2;
        %------------------------------------------------------------------
    end

    % run the test and time
    tic
    Q = DTQP_DEFECTS_convolution_integral_type0(A,B,in,opts);
    toc
end