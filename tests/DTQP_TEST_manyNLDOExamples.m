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
        opts.method.form = 'nonlinearprogram';
        opts.solver.function = 'ipfmincon';
        %------------------------------------------------------------------
        case 2
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        opts.dt.nt = 20; % number of time points
        opts.method.form = 'nonlinearprogram';
        opts.solver.function = 'ipfmincon';
        %------------------------------------------------------------------
        case 3
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        opts.solver.maxiters = 10;
        opts.solver.display = 'none';
        opts.method.form = 'qlin';
        opts.method.trustregionflag = false;
        opts.method.sqpflag = false;
        opts.method.delta = inf;
        opts.method.improveguess = false; % disabled
        opts.general.displevel = 1;
        %------------------------------------------------------------------
        case 4
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 100; % number of time points
        opts.method.form = 'nonlinearprogram';
        opts.solver.function = 'ipfmincon';
        opts.method.derivatives = 'complex'; % complex-step finite differencing
        %------------------------------------------------------------------
        case 5
        opts.dt.defects = 'TR';
        opts.dt.quadrature = 'CTR';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 20; % number of time points
        opts.method.form = 'nonlinearprogram';
        opts.solver.function = 'ipfmincon';
        opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
        opts.dt.meshr.tolerance = 1e-4;
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
BatchFermentorPenicillin(p,opts)
Brachistochrone(p,opts)
BrysonHo59(p,opts)
BrysonHo63(p,opts)
BrysonHo64(p,opts)
ChemicalReactor(p,opts)
ContainerCrane(p,opts)
FreeFlyingRobot(p,opts)
% HagerHouRao1(p,opts)  % fails with qlin (divide by 0)
HangGlider(p,opts)
HIVImmunology(p,opts)
Hypersensitive(p,opts)
Multiextremal(p,opts)
NeuenhofenKerriganX1(p,opts)
NeuenhofenKerriganX2(p,opts)
Nonlinear1D(p,opts)
OptimalConsumption(p,opts)
SecondOrderSingular(p,opts)
SimpleCoDesignTransfer(p,opts)
SimpleSASA(p,opts)
SimpleSuspensionSimultaneous(p,opts)
SpaceShuttleReentry(p,opts)
Train(p,opts)
TransferMaxRadius(p,opts)
TransferMinFuel(p,opts)
Tuberculosis(p,opts)
Tumor(p,opts)
TwoLinkRobot(p,opts)
Vanderpol(p,opts)

end