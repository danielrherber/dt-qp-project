%--------------------------------------------------------------------------
% LQRScalarTransfer_opts.m
% User options function for LQRScalarTransfer example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = LQRScalarTransfer_opts

% test number
num = 1;

switch num

case 1
    opts.Defectmethod = 'HS';
    opts.Quadmethod = 'CQHS';
    opts.NType = 'CGL';
    p.nt = 2000; % number of nodes

end

end