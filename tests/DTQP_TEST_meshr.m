%--------------------------------------------------------------------------
% DTQP_TEST_meshr.m
% Testing different mesh refinement schemes that are currently available
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:4;
% tests = [1,2];

% problem to use
problem = @BrysonHo166;
% problem = @BrysonDenham;

% go through the tests
for k = 1:length(tests)

    clear opts

    % options
    opts.general.plotflag = 0;
    opts.qp.disp = 'none';
    opts.dt.nt = 3; % number of nodes
    opts.dt.meshr.tolerance = 1e-3; % error tolerance

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        opts.dt.defects = 'HS';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.meshr.method = 'none';
        %------------------------------------------------------------------
        case 2
        opts.dt.defects = 'HS';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.meshr.method = 'richardson-doubling';
        %------------------------------------------------------------------
        case 3
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.meshr.method = 'richardson-doubling';
        %------------------------------------------------------------------
        case 4
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        opts.dt.meshr.method = 'richardson-doubling';
    end

    % run the test and time
    O = problem([],opts);

    % test analysis
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        disp(strcat("-----> ",string(O(1).value)," <-----"))
        %------------------------------------------------------------------
        case 2
        disp(strcat("-----> ",string(O(1).value)," <-----"))
        %------------------------------------------------------------------
        case 3
        disp(strcat("-----> ",string(O(1).value)," <-----"))
        %------------------------------------------------------------------
        case 4
        disp(strcat("-----> ",string(O(1).value)," <-----"))
    end
end