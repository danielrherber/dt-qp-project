%--------------------------------------------------------------------------
% BrysonHo116_plot.m
% Plot function for BrysonHo116 example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function BrysonHo116_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% plot state
solutionflag = true; %#ok<NASGU>
flag = 'plot-state'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
solutionflag = true; %#ok<NASGU>
flag = 'plot-control'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% error
figure('Color',wcolor); hold on

% plot error
flag = 'plot-error'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-error'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% mesh refinement
if isfield(in,'meshr')

    % step size plot
    DTQP_MESH_plotStepsize(in.meshr.T)

    % mesh error plot
    DTQP_MESH_plotError(in.meshr.T,in.meshr.errors,in.meshr.tol)

end

end