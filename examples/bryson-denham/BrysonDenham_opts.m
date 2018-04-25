%--------------------------------------------------------------------------
% BrysonDenham_opts.m
% User options function for BrysonDenham example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = BrysonDenham_opts

% test number
num = 1;

switch num

case 1
    % default parameters
    opts.plotflag = 1; % create the plots
    opts.saveflag = 0;
    opts.displevel = 2;
    opts.Defectmethod = 'HS';
    opts.Quadmethod = 'CQHS';
    opts.NType = 'ED';
    p.nt = 16; % number of nodes
    opts.reorder = 0;
    opts.solver = 'built-in';
    opts.tolerance = 1e-15;
    opts.maxiters = 200;
    opts.disp = 'iter';
    
case 2
    opts.Defectmethod = 'HS';
    opts.Quadmethod = 'CQHS';
    opts.NType = 'ED';
    p.nt = 4; % number of nodes
    
end

end