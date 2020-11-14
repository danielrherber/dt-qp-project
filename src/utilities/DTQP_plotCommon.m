%--------------------------------------------------------------------------
% DTQP_plotCommon.m
% Common plotting tasks
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
switch flag
    %----------------------------------------------------------------------
    case 'preliminary'
    % requires:
    if exist('closeflag','var')
        if closeflag
            close all
        end
    else
        close all
    end
    fontsize = 16;
    wcolor = [1 1 1]; % white color
    bcolor = [0 0 0]; % black color
    set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
    set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
    set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
    %----------------------------------------------------------------------
    case 'axis' % requires: bcolor, fontsize
    ha = gca; % get current axis handle
    ha.XAxis.Color = bcolor; % change the x axis color to black (not a dark gray)
    ha.YAxis.Color = bcolor; % change the y axis color to black (not a dark gray)
    ha.XAxis.FontSize = fontsize-6; % change x tick font size
    ha.YAxis.FontSize = fontsize-6; % change y tick font size
    ha.XAxis.Label.FontSize = fontsize; % change x label font size
    ha.YAxis.Label.FontSize = fontsize; % change y label font size
    ha.Box = 'on'; % box on
    ha.Layer = 'top'; % place the axes on top of the data
    ha.LineWidth = 1; % increase axis line width
    %----------------------------------------------------------------------
    case 'legend' % requires: bcolor, fontsize
    hl.FontSize = fontsize-16; % change legend font size
    hl.EdgeColor = bcolor; % change the legend border to black (not a dark gray)
    %----------------------------------------------------------------------
    case 'save' % requires: opts, pathplots, figname
    if opts.general.saveflag
        % save a pdf and png version
        filename = fullfile(pathplots,figname);
        str = ['export_fig ''',filename,''' -png -pdf'];
        eval(str)
    end
    %----------------------------------------------------------------------
    case 'plot-state' % requires: T, Y, sol, fontsize, solutionflag
    % line colors
    cArray = lines(size(Y,2));

    % plot state
    for i = 1:size(Y,2)
        plot(T,Y(:,i),'.','color',cArray(i,:),'markersize',12);
        if solutionflag
            plot(sol(2).T,sol(2).Y(:,i),'linewidth',2,'color',cArray(i,:));
        end
    end

    % axis
    xlabel('$t$ (time)','fontsize',fontsize)
    ylabel('$\xi$ (states)','fontsize',fontsize)

    % legend
    Lv = {};
    for i = 1:size(Y,2)
        Lv{end+1} = ['$\xi^{DT}_{',num2str(i),'}$'];
        if solutionflag
            Lv{end+1} = ['$\xi^*_{',num2str(i),'}$'];
        end
    end
    hL = legend(Lv);
    set(hL,'interpreter','latex','location','best',...
        'fontsize',fontsize-4)
    %----------------------------------------------------------------------
    case 'plot-control' % requires: T, U, sol, fontsize, solutionflag
    % line colors
    cArray = lines(size(U,2));

    % plot control
    for i = 1:size(U,2)
        plot(T,U(:,i),'.','color',cArray(i,:),'markersize',12);
        if solutionflag
            plot(sol(2).T,sol(2).U(:,i),'linewidth',2,'color',cArray(i,:));
        end
    end

    % axis
    xlabel('$t$ (time)','fontsize',fontsize)
    ylabel('$u$ (controls)','fontsize',fontsize)

    % legend
    Lv = {};
    for i = 1:size(U,2)
        Lv{end+1} = ['$u^{DT}_{',num2str(i),'}$'];
        if solutionflag
            Lv{end+1} = ['$u^*_{',num2str(i),'}$'];
        end
    end
    hL = legend(Lv);
    set(hL,'interpreter','latex','location','best',...
        'fontsize',fontsize-4)
    %----------------------------------------------------------------------
    case 'plot-error' % requires: T, Y, U, sol, fontsize, solutionflag
    % initialize
    Hv = []; Lv = {};

    % plot state error
    for i = 1:size(Y,2)
        ph1 = semilogy(T,abs(sol(1).Y(:,i)-Y(:,i))+1e-20,'linewidth',2);
        Hv = [Hv, ph1];
        Lv{i} = ['$\mathrm{log}_{10}\|\xi^*_{',num2str(i),'} - \xi^{DT}_{',num2str(i),'}\|$'];
    end

    % plot control error
    for i = 1:size(U,2)
        ph1 = semilogy(T,abs(sol(1).U(:,i)-U(:,i))+1e-20,'linewidth',2);
        Hv = [Hv, ph1];
        Lv{size(Y,2)+i} = ['$\mathrm{log}_{10}\|u^*_{',num2str(i),'} - u^{DT}_{',num2str(i),'}\|$'];
    end

    % axis
    ylim([1e-17 1e0])
    xlabel('$t$ (time)')
    ylabel('$\mathrm{log}_{10}$ error')
    ha = gca;
    ha.YScale = 'log';

    % legend
    hL = legend(Hv,Lv);
    set(hL,'location','best','fontsize',fontsize-4)

end