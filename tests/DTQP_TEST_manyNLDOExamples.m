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

    auxdata = [];
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
    problem(auxdata,opts)

    % test analysis

end

% problem structure
function problem(~,opts)

% run all examples
AlpRider([],opts)
BatchFermentorPenicillin([],opts)
Brachistochrone([],opts)
BrysonHo59([],opts)
BrysonHo63([],opts)
BrysonHo64([],opts)
ChemicalReactor([],opts)
ContainerCrane([],opts)
DynamicSoaring([],opts)
% EarthLaunch([],opts)
FreeFlyingRobot([],opts)
% HagerHouRao1([],opts)  % fails with qlin (divide by 0)
HangGlider([],opts)
HIVImmunology([],opts)
Hypersensitive([],opts)
% MineExtraction([],opts) % implementation does support symbolic derivatives
MountainCar([],opts)
Multiextremal([],opts)
NeuenhofenKerriganX1([],opts)
NeuenhofenKerriganX2([],opts)
Nonlinear1D([],opts)
OptimalConsumption([],opts)
SecondOrderSingular([],opts)
SimpleCoDesignTransfer([],opts)
SimpleSASA([],opts)
SimpleSuspensionSimultaneous([],opts)
SpaceShuttleReentry([],opts)
TimeEnergyInterceptor([],opts)
Train([],opts)
TransferMaxRadius([],opts)
TransferMinFuel([],opts)
TransferMinTime([],opts)
Tuberculosis([],opts)
% Tumor([],opts) % fails with qlin (need to investigate)
TwoLinkRobot([],opts)
Vanderpol([],opts)

end