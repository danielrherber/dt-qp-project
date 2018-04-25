%--------------------------------------------------------------------------
% Hager1_opts.m
% User options function for Hager1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = Hager1_opts

% test number
num = 1;

switch num

case 1
    opts.Defectmethod = 'PS';
    opts.Quadmethod = 'G';
    opts.NType = 'LGL';
    p.nt = 10;

case 2
    opts.Quadmethod = 'CQHS';
    opts.Defectmethod = 'HS';
    opts.NType = 'ED';
    p.nt = 10;

case 3
    opts.Quadmethod = 'CTR';
    opts.Defectmethod = 'TR';
    opts.NType = 'ED';
    p.nt = 10;

end

end