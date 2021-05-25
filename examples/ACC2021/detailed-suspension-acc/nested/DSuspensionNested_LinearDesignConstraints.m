%--------------------------------------------------------------------------
% DSuspensionNested_LinearDesignConstraints.m
% Linear outer-loop (plant design) constraints
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [A,b] = DSuspensionNested_LinearDesignConstraints(param)

% extract problem parameters
L0max = param.L0max; DoMax = param.DoMax;
dsc = param.dsc; dwt = param.dwt;
ld1 = param.ld1; ld2 = param.ld2;

% matrix
A = [4/DoMax,-1/DoMax,0,0,0,0,0; % spring manufacturing constraint
    -12/DoMax,1/DoMax,0,0,0,0,0; % spring manufacturing constraint
    1/DoMax,1/DoMax,0,0,0,0,0; % spring packaging constraint
    1/(2*(dsc+dwt)),-1/(2*(dsc+dwt)),0,0,0,1/(2*(dsc+dwt)),0; % suspension stop constraint
    0,0,0,0,0,0,2/(-ld1 - ld2 + L0max)]; % ensures enough space for damper

% vector
b = [0;0;1;-1;1];

end