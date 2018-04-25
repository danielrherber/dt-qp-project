%--------------------------------------------------------------------------
% DTQP2_opts.m
% User options function for DTQP2 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [p,opts] = DTQP2_opts

% test number
num = 1;

switch num

case 1
    opts = []; % defaults
    p.nt = 1000; % number of time points

end

end