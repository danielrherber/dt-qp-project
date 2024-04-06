%--------------------------------------------------------------------------
% TetherAssistedReentry_plot.m
% Plot function for TetherAssistedReentry example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function TetherAssistedReentry_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% plot state
solutionflag = false; %#ok<NASGU>
flag = 'plot-state'; DTQP_plotCommon %#ok<NASGU>

% modify lines
ha = gca;
ha.Children(1).LineStyle = '-';
ha.Children(1).LineWidth = 1.5;
ha.Children(2).LineStyle = '-';
ha.Children(2).LineWidth = 1.5;
ha.Children(3).LineStyle = '-';
ha.Children(3).LineWidth = 1.5;
ha.Children(4).LineStyle = '-';
ha.Children(4).LineWidth = 1.5;

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
solutionflag = false; %#ok<NASGU>
flag = 'plot-control'; DTQP_plotCommon %#ok<NASGU>

% modify lines
ha = gca;
ha.Children(1).LineStyle = '-';
ha.Children(1).LineWidth = 1.5;

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end