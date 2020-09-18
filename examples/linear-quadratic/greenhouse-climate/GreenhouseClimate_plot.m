%--------------------------------------------------------------------------
% GreenhouseClimate_plot.m
% Plot function for GreenhouseClimate example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function GreenhouseClimate_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract problem parameters structure
p = in.p;

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% plot state
solutionflag = false; %#ok<NASGU>
flag = 'plot-state'; DTQP_plotCommon %#ok<NASGU>

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

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% additional plots

% case-specific parameters
switch p.casenum
    %----------------------------------------------------------------------
    case {1,3}
    smax = 41;
    umax = 35;
    d1 = 0.5*(800*sin(4*pi*T/p.xtf - 0.65*pi) + abs(800*sin(4*pi*T/p.xtf - 0.65*pi)));
    us = 5;
    %----------------------------------------------------------------------
    case {2,4}
    smax = 30;
    umax = 30;
    d1 = 800*sin(4*pi*T/p.xtf - 0.65*pi);
    us = 1;
    %----------------------------------------------------------------------
    otherwise
        return % don't create the extra plots
end

%--- inputs and control
figure('Color',wcolor); hold on

% plot trajectories
s_light = 40;
plot(T,d1/s_light,'linewidth',2);
plot(T,15 + 10*sin(4*pi*T/p.xtf - 0.65*pi),'linewidth',2);
plot(T,U/us,'linewidth',2)

% labels
xlabel('Time [h]');
ylabel('Heat input, outside temperature, \& light');
legend(['Light$\times$',num2str(s_light),' [W]'],'Outside temp. [$^\circ$C]',...
    ['(Heat input)$/$',num2str(us),' [W]'],...
    'Location','best');

% configure axis
axis([0 p.xtf -1 umax]);
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-inputs-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%--- states and running cost
figure('Color',wcolor); hold on

% plot parameters
sf1 = 1200; sf3 = 60;

% plot trajectories
plot(T,[sf1*Y(:,1) Y(:,2) sf3*cumtrapz(T,p.xp4*U(:,1))],'linewidth',2);

% labels
xlabel('Time [h]'); ylabel('states');
legend('1200$\times$Dry weight [kg]','Greenhouse temp. [$^\circ$C]','$60 p_4 \int_{0}^{t_f}u\ dt$ [J]',...
    'Location','best');

% configure axis
axis([0 p.xtf -5 smax]);
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-states-cost'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end