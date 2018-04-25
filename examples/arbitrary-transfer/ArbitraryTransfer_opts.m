%--------------------------------------------------------------------------
% ArbitraryTransfer_opts.m
% User options function for ArbitraryTransfer example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = ArbitraryTransfer_opts

% test number
num = 1;

switch num

case 1
    opts.Defectmethod = 'HS';
    opts.Quadmethod = 'CQHS';
    opts.NType = 'ED'; 
    p.nt = 200; % number of time points
    
end

end