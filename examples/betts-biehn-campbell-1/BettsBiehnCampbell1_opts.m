%--------------------------------------------------------------------------
% BettsBiehnCampbell1_opts.m
% User options function for BettsBiehnCampbell1 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = BettsBiehnCampbell1_opts

% test number
num = 2;

switch num

case 1
    opts.Defectmethod = 'TR';
    opts.Quadmethod = 'CTR';
    opts.NType = 'ED';
    p.nt = 11; % number of nodes
    
case 2
    opts.Defectmethod = 'PS';
    opts.Quadmethod = 'G';
    opts.NType = 'LGL';
    p.nt = 11; % number of nodes
    
end

end