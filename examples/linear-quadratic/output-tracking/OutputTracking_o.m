%--------------------------------------------------------------------------
% OutputTracking_output.m
% Helper function for OutputTracking example
%--------------------------------------------------------------------------
% NOTE: can be an arbitrary function of time of the correct size and
% solution should automatically be created upon the next run
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [o,W] = OutputTracking_o(no)

W = linspace(0.1,1,no);
o = cell(no,1);
for k = 1:no
    o{k} = @(t) sin(W(k)*t);
end

end