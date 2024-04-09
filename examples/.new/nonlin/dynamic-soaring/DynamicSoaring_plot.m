%--------------------------------------------------------------------------
% DynamicSoaring_plot.m
% Plot function for DynamicSoaring example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DynamicSoaring_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

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

%% positions
figure('Color',wcolor); hold on

% extract
x = Y(:,1);
y = Y(:,2);
z = Y(:,3);

% plot
plot3(x,y,z,'-o','linewidth',2);

% change view
view(gca,[23 33.5]);

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% axis labels
xlabel('$x$','fontsize',fontsize_)
ylabel('$y$','fontsize',fontsize_)
zlabel('$z$','fontsize',fontsize_)

% z axis customization
ha = gca; % get current axis handle
ha.ZAxis.Color = bcolor; % change the x axis color to black (not a dark gray)
ha.ZAxis.FontSize = fontsize_-3; % change x tick font size
ha.ZAxis.Label.FontSize = fontsize_; % change z label font size
ha.LineWidth = 1;

% save
figname = 'figure-positions'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% path constraint
figure('Color',wcolor); hold on

% extract
auxdata = in.auxdata;

% compute path constraint
G = 0.5*auxdata.rho*auxdata.S/auxdata.m/auxdata.g*U(:,1).*Y(:,4).^2;

% plot
plot(T,G,'-','linewidth',2)
plot(T,auxdata.lmax*ones(size(T)),'--','linewidth',2)
plot(T,auxdata.lmin*ones(size(T)),'--','linewidth',2)

% axis limits
ylim([auxdata.lmin*1.1 auxdata.lmax*1.1])

% axis labels
xlabel('$t$ (time)','fontsize',fontsize_)
ylabel('$G$','fontsize',fontsize_)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-positions'; %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end