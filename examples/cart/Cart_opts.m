%--------------------------------------------------------------------------
% Cart_opts.m
% User options function for Cart example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = Cart_opts

% test number
num = 1;

switch num

case 1
    opts = [];
    p.nt = 100; % number of nodes

end

end