%--------------------------------------------------------------------------
% LQRstandard_plot.m
% Plot function for LQRstandard example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function LQRstandard_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

close all

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

fontsize = 16;

%% state
figure('Color',[1 1 1]);

% plot state
cArray = parula(size(Y,2));
for i = 1:size(Y,2)
    plot(T,Y(:,i),'.','color',cArray(i,:),'markersize',12); hold on
    plot(sol(2).T,sol(2).Y(:,i),'linewidth',2,'color',cArray(i,:)); hold on
end

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$\xi$','fontsize',fontsize)

% legend
Lv = {};
for i = 1:size(Y,2)
    Lv{end+1} = ['$\xi^{DT}_',num2str(i),'$'];
    Lv{end+1} = ['$\xi^*_',num2str(i),'$'];
end
hL = legend(Lv{1:2}); % only the first two entries
set(hL,'interpreter','latex','location','best','fontsize',fontsize-4,'box','on')

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-state'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

%% control
figure('Color',[1 1 1]);

% plot control
cArray = parula(size(U,2));
for i = 1:size(U,2)
    plot(T,U(:,i),'.','color',cArray(i,:),'markersize',12); hold on
    plot(sol(2).T,sol(2).U(:,i),'linewidth',2,'color',cArray(i,:)); hold on
end

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$u$','fontsize',fontsize)

% legend
Lv = {};
for i = 1:size(U,2)
    Lv{end+1} = ['$u^{DT}_',num2str(i),'$'];
    Lv{end+1} = ['$u^*_',num2str(i),'$'];
end
hL = legend(Lv{1:2});
set(hL,'interpreter','latex','location','best','fontsize',fontsize-4,'box','on')

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-control'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

%% error
figure('Color',[1 1 1]);

% initialize
Hv = []; Lv = {};

% plot state error
ph1 = semilogy(T,max(abs(sol(1).Y-Y),[],2)+1e-20,'linewidth',2); hold on
Lv{1} = ['$\mathrm{log}_{10}\max_k\|\xi^*_{k}(t) - \xi^{DT}_{k}(t) \|$'];

% plot control error
ph1 = semilogy(T,max(abs(sol(1).U-U),[],2)+1e-20,'linewidth',2); hold on
Lv{2} = ['$\mathrm{log}_{10}\max_k\|u^*_{k}(t) - u^{DT}_{k}(t) \|$'];

% axis
ylim([1e-17 1e0])
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$\mathrm{log}_{10}$ error','fontsize',fontsize)

% legend
hL = legend(Hv,Lv);
set(hL,'location','best','fontsize',fontsize-4,'box','on')

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-error'];
    str = ['export_fig ''',filename,''' -png -pdf'];
    eval(str)
end

end