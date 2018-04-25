%--------------------------------------------------------------------------
% HDAE_plot.m
% Plot function for HDAE example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function HDAE_plot(T,U,Y,P,F,p,opts,sol)

if opts.plotflag

close all

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

fontsize = 16;

%% state
figure('Color',[1 1 1]);

% plot state
cArray = flipud(parula(size(Y,2)));
for i = 1:size(Y,2)
    plot(T,Y(:,i),'.-','color',cArray(i,:),'markersize',12); hold on
%     plot(sol(2).T,sol(2).Y(:,i),'linewidth',2,'color',cArray(i,:)); hold on
end

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$\xi$','fontsize',fontsize)

% legend
% Lv = {};
% for i = 1:size(Y,2)
%     Lv{end+1} = ['$\xi^{DT}_{',num2str(i),'}$'];
% %     Lv{end+1} = ['$\xi^*_',num2str(i),'$'];
% end
% hL = legend(Lv{end}); % only the last entry
% set(hL,'interpreter','latex','location','best','fontsize',fontsize-4,'box','on')

% save
if opts.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-state'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

%% control
figure('Color',[1 1 1]);

% plot control
cArray = flipud(parula(size(U,2)));
for i = 1:size(U,2)
    plot(T,U(:,i),'.-','color',cArray(i,:),'markersize',12); hold on
%     plot(sol(2).T,sol(2).U(:,i),'linewidth',2,'color',cArray(i,:)); hold on
end

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$u$','fontsize',fontsize)

% legend
% Lv = {};
% for i = 1:size(U,2)
%     Lv{end+1} = ['$u^{DT}_{',num2str(i),'}$'];
% %     Lv{end+1} = ['$u^*_',num2str(i),'$'];
% end
% hL = legend(Lv{end}); % only the last entry
% set(hL,'interpreter','latex','location','best','fontsize',fontsize-4,'box','on')

% save
if opts.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-control'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

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