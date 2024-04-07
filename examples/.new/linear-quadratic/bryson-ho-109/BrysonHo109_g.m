%--------------------------------------------------------------------------
% BrysonHo109_g.m
% Helper function for BrysonHo109 example
%--------------------------------------------------------------------------
% NOTE: can be an arbitrary function of time and solution should
% automatically be created upon the next run
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function g = BrysonHo109_g(t)
g = t.*cos(20*pi*t) - 1/4;
end