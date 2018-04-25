%--------------------------------------------------------------------------
% HDAE_opts.m
% User options function for HDAE example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = HDAE_opts

% test number
num = 1;

switch num

case 1
    opts.Defectmethod = 'HS';
    opts.Quadmethod = 'CQHS';
    opts.NType = 'ED';
    p.nt = 100;

case 2
    opts.Defectmethod = 'PS';
    opts.Quadmethod = 'G';
    opts.NType = 'LGL';
    p.nt = 100; % number of nodes

end

end