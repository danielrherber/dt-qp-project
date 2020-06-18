%--------------------------------------------------------------------------
% DTQP_TEST_ipfmincon_symb.m
% Tests for DTQP_ipfmincon_symb.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
clc; clear; close all

tests = 1:9;
% tests = 1;

% go through the tests
for k = 1:length(tests)

    clear f o form
    
    % options
    D2flag = true;
    linflag = true;
    opts.general.displevel = 0;

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        % linearization of the quadratic
        f = 'y1*(1-y1)';
        o.ny = 1; % number of states
        o.nu = 0; % number of controls
        %------------------------------------------------------------------
        case 2
        % quadraticization of the cubic ( additional symbolic parameters)
        f = 'a*b*y1*(1-y1)^2';
        o.ny = 1; % number of states
        o.paramstr = 'a b';
        %------------------------------------------------------------------
        case 3
        % linearization of the Van der Pol oscillator equations of motion
        f = '[y2;-y1 + y2 - y1^2*y2 + u1]';
        o.ny = 2; % number of states
        o.nu = 1; % number of controls
        %------------------------------------------------------------------
        case 4
        % linearization of the Hyper-Sensitive equations of motion
        f = '-y1^3 + u1';
        o.ny = 1; % number of states
        o.nu = 1; % number of controls
        %------------------------------------------------------------------
        case 5
        %
        f = '[alp*y1*(1+sqrt(1-y1))/2]';
        o.ny = 3; % number of states
        o.nu = 1; % number of controls
        o.np = 0; % number of parameters
        o.output = 2; % interpolated version
        o.paramstr = 'alp';
        %------------------------------------------------------------------
        case 6
        %
        f = '[sin(t)*y1^3+y2^2+y1*y2+y1+y2]';
        o.ny = 2; % number of states
        o.nu = 0; % number of controls
        o.np = 0; % number of parameters
        o.output = 2; % interpolated version
        %------------------------------------------------------------------
        case 7
        f = '[y1^2*y2]';
        o.ny = 2; % number of states
        o.nu = 0; % number of controls
        o.np = 0; % number of parameters
        o.output = 2; % interpolated version
        %------------------------------------------------------------------
        case 8
        %
        f = '[-zeta.*y1.*log(y1./y2); y2.*(b-(Mew+(D*(y1.^(2/3)))+G.*u1))]';
        o.ny = 2; % number of states
        o.nu = 1; % number of controls
        o.np = 0; % number of parameters
        o.paramstr = 'zeta b Mew D G';
        %------------------------------------------------------------------
        case 9
        %
        f = '[zeta*sin(t)*u1;cos(t)*t^2*y1;p1+u1+y2+t]';
        o.ny = 3; % number of states
        o.nu = 1; % number of controls
        o.np = 1; % number of parameters
        o.paramstr = 'zeta';
    end

    % run the test and time
    tic
    E = DTQP_ipfmincon_symb(f,o,linflag,opts);
    toc

    % test analysis
    disp(E); disp(' ');

end