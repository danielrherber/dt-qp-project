%--------------------------------------------------------------------------
% TimeEnergyInterceptor_plot.m
% Plot function for TimeEnergyInterceptor example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function TimeEnergyInterceptor_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract final time;
tf = P(1);
T = tf*T;

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
figname = 'figure-state'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
solutionflag = false; %#ok<NASGU>
flag = 'plot-control'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% axis 1 vs. axis 2
if in.auxdata.formulation == 1
figure('Color',wcolor); hold on

% extract radius
r = in.auxdata.r;

% plot trajectories
plot(Y(:,2),Y(:,3),'b.-','LineWidth',1,'markersize',12,...
    'Color',[165, 61, 151]/255,'DisplayName',"Target")
plot(Y(:,5),Y(:,6),'k.-','LineWidth',1,'markersize',12,'DisplayName',"Missile")

% plot final lethal radius
rectangle('Position',[Y(end,2)-r Y(end,3)-r 2*r 2*r],...
    'Curvature',[1,1],'LineWidth',1,'EdgeColor','b');

% legend
hl = legend;
flag = 'legend'; DTQP_plotCommon %#ok<NASGU>

% axis labels
xlabel('axis 1 [m]')
ylabel('axis 2 [m]')

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-axis1-axis2'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end

end