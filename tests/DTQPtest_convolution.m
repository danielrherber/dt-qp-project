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
close all
clear
clc

testnum = 1;

switch testnum
%--------------------------------------------------------------------------
case 1 % test 1
    opts.NType = 'ED';
    p.nt = 23;
    p.t = linspace(0,1,p.nt);
    p.h = diff(p.t);

    % constant matrices
    A = rand(3);
    B = [1 1;1 2;1 3];
    
    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
case 2 % test 2
    opts.NType = 'ED';
    p.nt = 10000;
    p.t = linspace(0,1,p.nt);
    p.h = diff(p.t);

    % singular A
    A = [0 1; 0 0];
    B = [1;1];

    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
case 3 % test 3
    opts.NType = 'LGL';
    p.nt = 100;
    p.t = DTQP_nodes_LGL(p.nt-1);
    p.h = diff(p.t);

    % time-varying matrix
    A = magic(2);
    B{1,1} = 0;
    B{2,1} = @(t) exp(-t);
    
    % time
    tic
    Q = DTQP_convolution(A,B,p,opts);
    toc
%--------------------------------------------------------------------------
end