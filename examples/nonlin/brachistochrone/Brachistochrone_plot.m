%--------------------------------------------------------------------------
% Brachistochrone_plot.m
% Plot function for Brachistochrone example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function Brachistochrone_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract parameter structure
p = in.p;

% solution available for all cases
solutionflag = true;

% unscale time
T = T*P(1);

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% plot state
flag = 'plot-state'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>
xlim([T(1),max(sol(1).F,T(end))])

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
flag = 'plot-control'; DTQP_plotCommon %#ok<NASGU>

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>
xlim([T(1),max(sol(1).F,T(end))])

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% error
if solutionflag
    figure('Color',wcolor); hold on

    % plot error
    flag = 'plot-error'; DTQP_plotCommon %#ok<NASGU>

    % configure axis
    flag = 'axis'; DTQP_plotCommon %#ok<NASGU>
    xlim([T(1),max(sol(1).F,T(end))])

    % save
    figname = 'figure-error'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
    flag = 'save'; DTQP_plotCommon %#ok<NASGU>
end

%% x-y
figure('Color',wcolor); hold on

% initialize legend
Lv = {};

% line colors
cArray = lines(size(Y,2));

% hard-coded path constraint
if p.casenum == 4
    plot([min(Y(:,1)),max(Y(:,1))],-[min(Y(:,1))/2 + 0.1, max(Y(:,1))/2 + 0.1],...
        '-','linewidth',2,'color',cArray(2,:));
    Lv{end+1} = ['path constraint'];
end

% plot state
plot(Y(:,1),-Y(:,2),'.','color',cArray(1,:),'markersize',12);
if solutionflag
    plot(sol(2).Y(:,1),-sol(2).Y(:,2),'linewidth',2,'color',cArray(1,:));
end

% axis
xlabel('$x$','fontsize',fontsize)
ylabel('$y$','fontsize',fontsize)

% legend
Lv{end+1} = ['$(x,y)^{DT}$'];
if solutionflag
    Lv{end+1} = ['$(x,y)^*$'];
end
hL = legend(Lv);
set(hL,'interpreter','latex','location','best',...
    'fontsize',fontsize-4)

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-xy'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

end