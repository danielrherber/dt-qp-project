%--------------------------------------------------------------------------
% LQRInhomogeneous_d.m
% Helper function for LQRInhomogeneous example
%--------------------------------------------------------------------------
% NOTE: can be an arbitrary function of time of the correct size and
% solution should automatically be created upon the next run
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function d = LQRInhomogeneous_d(ny)

d = cell(ny,1);
d{1} = @(t) 10*sin(3*t);
d{3} = @(t) 2*sin(2*t);
d{4} = @(t) 4*sin(6*t);
d{4} = @(t) -3*sin(1*t);
d{5} = @(t) -1*sin(2*t);
d{6} = @(t) 10*sin(3*t);
d{7} = @(t) 2*sin(2*t);
d{8} = @(t) 4*sin(6*t);
d{9} = @(t) -3*sin(1*t);
d{10} = @(t) -1*sin(2*t);

end