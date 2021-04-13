%--------------------------------------------------------------------------
% DTQP_DEFECTS_expm.m
% Matrix exponential computed only for the unique step sizes
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [At,At2,At3] = DTQP_DEFECTS_expm(A,in,~)

% check A type
if ~isa(A,'double')
    error('A should not be a function or time varying with the zero-order hold (ZOH) or first-order hold (FOH) methods')
end

% extract
h = in.h;

% size of square A matrix
na = size(A,2);

% tolerance for uniquetol (NOTE: potentially expose)
tol = 1e-12; % default

% find unique step size values
[h_unique,~,IC_unique] = uniquetol(h,tol);

% number of unique step sizes
nh = length(h_unique);

% initialize
At_unique = zeros(nh,na,na);

% go through each unique step size
for k = 1:nh

    % compute matrix exponential
    At_unique(k,:,:) = expm(A*h_unique(k));

end

% assign
At = At_unique(IC_unique,:,:);

end