%--------------------------------------------------------------------------
% DSuspensionNested_InnerLoopRamp.m
% Inner-loop problem for the ramp load case
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [F,varargout] = DSuspensionNested_InnerLoopRamp(x,auxdata)

% common inner-loop problem setup
[setup,Bz] = DSuspensionNested_InnerLoopCommonSetup(x,auxdata);

% dynamic constraints
setup = DSuspensionNested_RampStateConstraints(x,auxdata,setup);

% combine
setup.t0 = 0; setup.tf = auxdata.tf; setup.auxdata = auxdata;

% ramp disturbance
ns = size(Bz,1);
d = cell(ns,1);
for idx = 1:ns
    d{idx,1} = @(t,p) Bz(idx)*p.ramp_in;
end
setup.d = d;

% DTQP options
opts.general.displevel = 0;
opts.general.plotflag = 0;
opts.solver.tolerance = auxdata.InnerLoopTolerance;
opts.dt.defects = 'TR';
opts.dt.quadrature = 'CTR';
opts.dt.mesh = 'ED';
opts.dt.nt = auxdata.nt;

% form and solve problem
[T,U,Y,P,F,in,opts] = DTQP_solve(setup,opts);

% optional additional outputs
if nargout > 1
    out.T = T; out.U = U; out.Y = Y; out.P = P; out.F = F;
    out.in = in; out.opts = opts;
    varargout{1} = out;
end

end