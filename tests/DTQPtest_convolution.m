%--------------------------------------------------------------------------
% DTQPtest_convolution.m
% Test the computation of the convolution integral
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

testnum = 1;

switch testnum
%--------------------------------------------------------------------------
case 1 % nonsingular A, ED mesh
    opts.dt.mesh = 'ED';
    opts.dt.nt = 23;
    p.t = linspace(0,1,opts.dt.nt);
    p.h = diff(p.t);

    % constant matrices
    A = rand(3);
    B = [1 1;1 2;1 3];
    
    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
case 2 % nonsingular A, USER mesh
    opts.dt.mesh = 'USER';
    p.t = [0,0.1,0.5,0.7,1];
    p.h = diff(p.t);
    opts.dt.nt = length(p.t);

    % constant matrices
    A = rand(3);
    B = [1 1;1 2;1 3];
    
    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
    
%--------------------------------------------------------------------------
case 3 % singular A, ED mesh
    opts.dt.mesh = 'ED';
    opts.dt.nt = 10000;
    p.t = linspace(0,1,opts.dt.nt);
    p.h = diff(p.t);

    % singular A
    A = [0 1; 0 0];
    B = [1;1];

    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
case 4 % singular A, USER mesh
    opts.dt.mesh = 'USER';
    p.t = [0,0.1,0.5,0.7,1];
    p.h = diff(p.t);
    opts.dt.nt = length(p.t);

    % singular A
    A = [0 1; 0 0];
    B = [1;1];

    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
case 5 % time-varying B
    opts.dt.mesh = 'LGL';
    opts.dt.nt = 100;
    p.t = DTQP_nodes_LGL(opts.dt.nt-1);
    p.h = diff(p.t);

    % time-varying matrix
    A = magic(2);
    B{1,1} = 0;
    B{2,1} = @(t) exp(-t);
    B{1,2} = @(t) sin(t);
    B{2,2} = 1;

    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
end