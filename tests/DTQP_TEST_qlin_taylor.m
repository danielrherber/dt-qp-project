%--------------------------------------------------------------------------
% DTQP_TEST_qlin_taylor.m
% A collection of test cases for DTQP_qlin_taylor.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

% go through each test
for k = 1:9
    tic
    disp(k)
    E = taylor_reduction_testcases(k);
    toc
    disp(E);
end

function E = taylor_reduction_testcases(testnum)
switch testnum
    %----------------------------------------------------------------------
    case 1
        clear
        % linearization of the quadratic
        f = 'y1*(1-y1)';
        o.ny = 1; % number of states
        o.nu = 0; % number of inputs
        form = 3;
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 2
        clear
        % quadraticization of the cubic (optional additional symbolic variables)
        f = 'a*b*y1*(1-y1)^2';
        o.ny = 1; % number of states
        o.param = 'a b';
        form = 4;
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 3
        clear
        % linearization of the Van der Pol oscillator equations of motion
        f = '[y2;-y1 + y2 - y1^2*y2 + u1]';
        o.ny = 2; % number of states
        o.nu = 1; % number of inputs
        form = 3;
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 4
        % linearization of the Hyper-Sensitive equations of motion
        f = '-y1^3 + u1';
        o.ny = 1; % number of states
        o.nu = 1; % number of inputs
        form = 3;
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 5
        clear
        % create
        f = '[alp*y1*(1+sqrt(1-y1))/2]';
        form = 4;
        o.ny = 3; % number of states
        o.nu = 1; % number of inputs
        o.np = 0; % number of parameters
        o.output = 2; % interpolated version
        o.param = 'alp';
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 6
        clear
        % create
        f = '[sin(t)*y1^3+y2^2+y1*y2+y1+y2]';
        form = 4;
        o.ny = 2; % number of states
        o.nu = 0; % number of inputs
        o.np = 0; % number of parameters
        o.output = 2; % interpolated version
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 7
        clear
        % create
        f = '[y1^2*y2]';
        form = 3;
        o.ny = 2; % number of states
        o.nu = 0; % number of inputs
        o.np = 0; % number of parameters
        o.output = 2; % interpolated version
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 8
        clear
        % create
        f = '[-zeta.*y1.*log(y1./y2); y2.*(b-(Mew+(D*(y1.^(2/3)))+G.*u1)); u1]';
        form = 3;
        o.ny = 2; % number of states
        o.nu = 1; % number of inputs
        o.np = 0; % number of parameters
        o.param = 'zeta b Mew D G';
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
    case 9
        clear
        % create
        f = '[zeta*sin(t)*u1;cos(t)*t^2*y1;p1+u1+y2]';
        form = 3;
        o.ny = 2; % number of states
        o.nu = 1; % number of inputs
        o.np = 1; % number of parameters
        o.param = 'zeta';
        E = DTQP_qlin_taylor(f,form,o);
    %----------------------------------------------------------------------
end
end