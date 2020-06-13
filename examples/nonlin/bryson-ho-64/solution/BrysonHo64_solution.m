%--------------------------------------------------------------------------
% BrysonHo64_solution.m
% Create solution for BrysonHo64 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [U,Y,F] = BrysonHo64_solution(a,l,T)

% find l/H
lH = fzero(@(x) a/l*x - cosh(x),0,optimset('Display','iter'));

% compute H
H = l/lH;

% calculate two candidate objectives (area)
F1 = 2*pi*a^2*(tanh(lH) + lH*sech(lH)^2);
F2 = 2*pi*a^2;

% determine which objective if smaller and calculate state and control
if F1 < F2
    F = F1;
    Y = H*cosh(T/H);
    U = sinh(T/H);
else
    F = F2;
    Y = repelem(a,length(T),1);
    U = zeros(size(T));
end

end