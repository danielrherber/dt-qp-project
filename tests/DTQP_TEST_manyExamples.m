%--------------------------------------------------------------------------
% DTQP_TEST_manyExamples.m
% Run many of the examples
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:4;
% tests = 1;

% go through each test
for k = 1:length(tests)

    p = [];
    opts.general.plotflag = 0;

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        opts.dt.defects = 'HS';
        opts.dt.quadrature = 'CQHS';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        %------------------------------------------------------------------
        case 2
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        opts.dt.nt = 20; % number of time points
        %------------------------------------------------------------------
        case 3
        opts.dt.defects = 'EF';
        opts.dt.quadrature = 'CEF';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        %------------------------------------------------------------------
        case 4
        opts.dt.defects = 'ZO';
        opts.dt.quadrature = 'CEF';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 20; % number of time points
    end

    % run the test and time
    problem(p,opts)

    % test analysis

end

% problem structure
function problem(p,opts)

% run all examples
AndersonMoore64(p,opts)
BettsBiehnCampbell1(p,opts)
BrysonDenham(p,opts)
BrysonHo109(p,opts)
BrysonHo116(p,opts)
BrysonHo153(p,opts)
BrysonHo154(p,opts)
BrysonHo156(p,opts)
BrysonHo166(p,opts)
BrysonHo248(p,opts)
Cart(p,opts)
DTQP1(p,opts)
DTQP2(p,opts)
DTQP3(p,opts)
Hager1(p,opts)
HDAE(p,opts)
JadduShimemura(p,opts)
LinearPendulum(p,opts)
LQRInhomogeneous(p,opts)
LQRScalar(p,opts)
LQRScalarTransfer(p,opts)
% LQRstandard(p,opts)
MinimumEnergyTransfer(p,opts)
MultiphaseParameter(p,opts)
% OutputTracking(p,opts)
TavallaeiTousi1(p,opts)
TurnerChunJuang1(p,opts)

end