%--------------------------------------------------------------------------
% DTQP1_opts.m
% User options function for DTQP1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = DTQP1_opts

% test number
num = 1;

switch num

case 1
    opts.plotflag = 1; % create the plots
    opts.saveflag = 1;
    opts.displevel = 2;
    opts.Defectmethod = 'HS';
    opts.Quadmethod = 'CQHS';
    opts.NType = 'ED';
    p.nt = 5000; % number of nodes
    opts.reorder = 0;
    opts.solver = 'built-in';
    opts.tolerance = 1e-15;
    opts.maxiters = 100;
    opts.disp = 'iter';

end

end