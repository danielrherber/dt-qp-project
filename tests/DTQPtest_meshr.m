%--------------------------------------------------------------------------
% DTQPtest_meshr.m
% Testing different mesh refinement schemes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
close all; clear; clc

problem = @BrysonHo166;
% problem = @BrysonDenham;

testnum = 1;

switch testnum
    
%--------------------------------------------------------------------------
    case 1
        opts.general.plotflag = 0;
        opts.dt.defects = 'HS';
        opts.dt.quadrature = 'CQHS';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 2; % number of nodes
        opts.dt.meshr.method = 'richardson-doubling';

        O = problem([],opts);

%--------------------------------------------------------------------------
    case 2
        opts.general.plotflag = 0;
        opts.dt.defects = 'PS';
        opts.dt.quadrature = 'G';
        opts.dt.mesh = 'LGL';
        % opts.dt.nt = 7; % number of nodes
        opts.dt.meshr.method = 'richardson-doubling';

        O = problem([],opts);

%--------------------------------------------------------------------------
    case 3
        opts.general.plotflag = 0;
        opts.dt.defects = 'EF';
        opts.dt.quadrature = 'CEF';
        opts.dt.mesh = 'ED';
        % opts.dt.nt = 7; % number of nodes
        opts.dt.meshr.method = 'richardson-doubling';

        O = problem([],opts);

%--------------------------------------------------------------------------
    case 4
        opts.dt.defects = 'HS';
        opts.dt.quadrature = 'CQHS';
        opts.dt.mesh = 'ED';
        opts.dt.nt = 20; % number of nodes
        opts.dt.meshr.method = 'none';

        O = problem([],opts);
end