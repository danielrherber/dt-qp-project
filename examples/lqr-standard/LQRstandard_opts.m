%--------------------------------------------------------------------------
% LQRstandard_opts.m
% User options function for LQRstandard example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = LQRstandard_opts

% test number
num = 1;

switch num

case 1
    % default parameters
    opts.plotflag = 1; % create the plots
    opts.saveflag = 0;
    opts.displevel = 2;
    opts.Defectmethod = 'PS';
    opts.Quadmethod = 'G';
    opts.NType = 'LGL';
    opts.reorder = 0;
    opts.solver = 'built-in';
    opts.tolerance = 1e-15;
    opts.maxiters = 200;
    opts.disp = 'iter';
    p.nt = 50; % number of nodes

end

end