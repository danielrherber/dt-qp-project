%--------------------------------------------------------------------------
% DTQP_TEST_manyNLDOExamples.m
% Run many of the NLDO examples
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

tests = 1:4;
% tests = 2;

% go through each test
for k = 1:length(tests)

    p = [];
    clear opts
    opts.general.plotflag = 0;

    % test setup
    switch tests(k)
        %------------------------------------------------------------------
        case 1
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        opts.qlin.method = 'ipfmincon';
        opts.qp.solver = 'ipfmincon';
        %------------------------------------------------------------------
        case 2
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        opts.dt.nt = 20; % number of time points
        opts.qlin.method = 'ipfmincon';
        opts.qp.solver = 'ipfmincon';
        %------------------------------------------------------------------
        case 3
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        opts.qp.maxiters = 10;
        opts.qp.disp = 'none';
        opts.qlin.method = 'qlin';
        opts.qlin.trustregionflag = false;
        opts.qlin.sqpflag = false;
        opts.qlin.delta = inf;
        opts.qlin.improveX0flag = false; % disabled
        opts.general.displevel = 1;
        %------------------------------------------------------------------
        case 4
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        opts.qlin.method = 'ipfmincon';
        opts.qp.solver = 'ipfmincon';
        opts.qlin.derivativemethod = 'complex'; % complex-step finite differencing
        %------------------------------------------------------------------
    end

    % run the test and time
    problem(p,opts)

    % test analysis

end

% problem structure
function problem(p,opts)

% run all examples
AlpRider(p,opts)
Brachistochrone(p,opts)
BrysonHo59(p,opts)
BrysonHo63(p,opts)
BrysonHo64(p,opts)
ChemicalReactor(p,opts)
ContainerCrane(p,opts)
FreeFlyingRobot(p,opts)
HangGlider(p,opts)
HIVImmunology(p,opts)
Hypersensitive(p,opts)
OptimalConsumption(p,opts)
SimpleCoDesignTransfer(p,opts)
Train(p,opts)
TransferMaxRadius(p,opts)
% TransferMinFuel(p,opts) % fails with qlin?
Tuberculosis(p,opts)
Tumor(p,opts)
Vanderpol(p,opts)

end