%--------------------------------------------------------------------------
% HDAE_plot.m
% Plot function for HDAE example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function HDAE_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract parameter structure
p = in.p;

% preliminary plot options
flag = 'preliminary'; DTQP_plotCommon %#ok<NASGU>

%% state
figure('Color',wcolor); hold on

% line colors
cArray = flipud(parula(size(Y,2)));

% plot state
for i = 1:size(Y,2)
    plot(T,Y(:,i),'.-','color',cArray(i,:),'markersize',12);
end

% axis
xlabel('$t$ (s)')
ylabel('$\xi$')

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-state'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% control
figure('Color',wcolor); hold on

% plot control
cArray = flipud(parula(size(U,2)));
for i = 1:size(U,2)
    plot(T,U(:,i),'.-','color',cArray(i,:),'markersize',12);
end

% axis
xlabel('$t$ (s)')
ylabel('$u$')

% configure axis
flag = 'axis'; DTQP_plotCommon %#ok<NASGU>

% save
figname = 'figure-control'; pathplots = msavename(mfilename('fullpath'),'plots'); %#ok<NASGU>
flag = 'save'; DTQP_plotCommon %#ok<NASGU>

%% states and controls (surf)
hf = figure('Color',[1 1 1]);
ha = axes('Parent',hf);
hold(ha,'on');

% combine
Z = [U(:,1),Y,U(:,end)];
TT = repmat(T,1,p.n+1);
x = (0:p.n)*(pi/p.n);
X = repmat(x,length(T),1);

% plot the data
surf(TT,X,Z,'Parent',ha);

% change view
view(ha,[-10.7 34]);

% change colormap
colormap(white)

% axis
ha.TickDir = 'in';
ha.XMinorTick = 'on';
ha.YMinorTick = 'on';
ha.ZMinorTick = 'on';
ha.XLim = [0 5];
ha.YLim = [0 4];
ha.ZLim = [-0.2 0.6];
ha.XLabel.String = 'time';
ha.YLabel.String = 'distance';
ha.ZLabel.String = 'temperature';
ha.YLabel.Rotation = -66;

end