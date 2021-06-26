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

tests = 1:5;
% tests = 1;

% go through each test
for k = 1:length(tests)

    auxdata = [];
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
        %------------------------------------------------------------------
        case 5
        opts.dt.defects = 'FO';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 20; % number of time points
    end

    % run the test and time
    problem(auxdata,opts)

    % test analysis

end

% problem structure
function problem(~,opts)

% run all examples
AndersonMoore64([],opts)
BettsBiehnCampbell1([],opts)
BrysonDenham([],opts)
BrysonHo109([],opts)
BrysonHo116([],opts)
BrysonHo153([],opts)
BrysonHo154([],opts)
BrysonHo156([],opts)
BrysonHo166([],opts)
BrysonHo248([],opts)
Cart([],opts)
DTQP1([],opts)
DTQP2([],opts)
DTQP3([],opts)
GasAbsorber([],opts)
% GreenhouseClimate([],opts) % fails (on purpose) with ZO because A(t)
Hager1([],opts)
HDAE([],opts)
JadduShimemura([],opts)
LinearPendulum([],opts)
LQRInhomogeneous([],opts)
LQRScalar([],opts)
LQRScalarTransfer([],opts)
LQRstandard([],opts)
MinimumEnergyTransfer([],opts)
MultiphaseParameter([],opts)
Nagurka([],opts)
OutputTracking([],opts)
TavallaeiTousi1([],opts)
TurnerChunJuang1([],opts)

end