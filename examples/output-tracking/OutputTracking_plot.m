%--------------------------------------------------------------------------
% OutputTracking_plot.m
% Plot function for OutputTracking example
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary Contributor: Daniel R. Herber, Graduate Student, University of 
% Illinois at Urbana-Champaign
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function OutputTracking_plot(T,U,Y,P,F,in,opts,sol)

if opts.general.plotflag

% extract parameter structure
p = in.p;

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
    plot(T,Y(:,i),'.','linewidth',2,'color',cArray(i,:),'markersize',12); hold on
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
    plot(T,U(:,i),'.','linewidth',2,'color',cArray(i,:),'markersize',12); hold on
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

%% output
figure('Color',[1 1 1]);

% plot state
cArray = parula(size(p.C,1));
for i = 1:size(p.C,1)
    plot(T, p.o{i}(T),'--','linewidth',2,'color',cArray(i,:)); hold on
    o(:,i) = sum(bsxfun(@times,p.C(i,:),Y),2);
    plot(T,o(:,i),'.','linewidth',2,'color',cArray(i,:),'markersize',12); hold on
    o_actual(:,i) = sum(bsxfun(@times,p.C(i,:),sol(2).Y),2);
    plot(sol(2).T,o_actual(:,i),'-','linewidth',2,'color',cArray(i,:),'markersize',12); hold on
end

% axis
xlabel('$t$ (s)','fontsize',fontsize)
ylabel('$o$ (output)','fontsize',fontsize)

% legend
Lv = {};
for i = 1:size(Y,2)
    Lv{end+1} = ['$y^*_',num2str(i),'$'];
    Lv{end+1} = ['$y^{DT}_',num2str(i),'$'];
end
hL = legend(Lv{1:2}); % only the first two entries
set(hL,'interpreter','latex','location','best','fontsize',fontsize-4,'box','on')

% save
if opts.general.saveflag
    path = msavename(mfilename('fullpath'),'plots');
    filename = [path,'figure-output'];
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