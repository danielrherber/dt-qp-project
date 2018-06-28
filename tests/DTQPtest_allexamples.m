%--------------------------------------------------------------------------
% DTQPtest_allexamples.m
% Run all the examples
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

p = [];

%% options 
opts.dt.defects = 'HS';
opts.dt.quadrature = 'CQHS';
opts.dt.mesh = 'ED'; 
opts.dt.nt = 1000; % number of time points
opts.general.plotflag = 0;

%% 
AndersonMoore64(p,opts)
ArbitraryTransfer(p,opts)
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
LQRInhomogeneous(p,opts)
LQRScalar(p,opts)
LQRScalarTransfer(p,opts)
LQRstandard(p,opts)
MultiphaseParameter(p,opts)
OutputTracking(p,opts)
