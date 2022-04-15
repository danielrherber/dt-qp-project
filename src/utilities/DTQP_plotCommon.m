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
    fontsize_ = 16;
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
    ha.XAxis.FontSize = fontsize_-3; % change x tick font size
    ha.YAxis.FontSize = fontsize_-3; % change y tick font size
    ha.XAxis.Label.FontSize = fontsize_; % change x label font size
    ha.YAxis.Label.FontSize = fontsize_; % change y label font size
    ha.Box = 'on'; % box on
    ha.Layer = 'top'; % place the axes on top of the data
    ha.LineWidth = 1; % increase axis line width
    %----------------------------------------------------------------------
    case 'legend' % requires: bcolor, fontsize
    hl.FontSize = fontsize_-3; % change legend font size
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

    % plot each state
    for i = 1:size(Y,2)

        % plot DT solution
        plot(T,Y(:,i),'.','markersize',12,'color',cArray(i,:),...
            'DisplayName',['$x^{DT}_{',num2str(i),'}$']);

        % plot alternative solution
        if solutionflag
            plot(sol(2).T,sol(2).Y(:,i),'linewidth',2,'color',cArray(i,:),...
                'DisplayName',['$x^*_{',num2str(i),'}$']);
        end

    end

    % axis
    xlabel('$t$ (time)','fontsize',fontsize_)
    ylabel('$x$ (states)','fontsize',fontsize_)

    % legend
    hL_opts = {'Interpreter','latex','FontSize',fontsize_-4,...
        'Location','best','EdgeColor',bcolor};
    legend(hL_opts{:});

    %----------------------------------------------------------------------
    case 'plot-control' % requires: T, U, sol, fontsize, solutionflag

    % line colors
    cArray = lines(size(U,2));

    % plot each control
    for i = 1:size(U,2)

        % plot DT solution
        plot(T,U(:,i),'.','markersize',12,'color',cArray(i,:),...
            'DisplayName',['$u^{DT}_{',num2str(i),'}$']);

        % plot alternative solution
        if solutionflag
            plot(sol(2).T,sol(2).U(:,i),'linewidth',2,'color',cArray(i,:),...
                'DisplayName',['$u^*_{',num2str(i),'}$']);
        end
    end

    % axis
    xlabel('$t$ (time)','fontsize',fontsize_)
    ylabel('$u$ (controls)','fontsize',fontsize_)

    % legend
    hL_opts = {'Interpreter','latex','FontSize',fontsize_-4,...
        'Location','best','EdgeColor',bcolor};
    legend(hL_opts{:});

    %----------------------------------------------------------------------
    case 'plot-error' % requires: T, Y, U, sol, fontsize, solutionflag

    % plot state error
    for i = 1:size(Y,2)
        semilogy(T,abs(sol(1).Y(:,i)-Y(:,i))+1e-20,'linewidth',2,...
            'DisplayName',['$\mathrm{log}_{10}\|x^*_{',num2str(i),'} - x^{DT}_{',num2str(i),'}\|$']);
    end

    % plot control error
    for i = 1:size(U,2)
        semilogy(T,abs(sol(1).U(:,i)-U(:,i))+1e-20,'linewidth',2,...
            'DisplayName',['$\mathrm{log}_{10}\|u^*_{',num2str(i),'} - u^{DT}_{',num2str(i),'}\|$']);
    end

    % axis
    ylim([1e-17 1e0])
    xlabel('$t$ (time)')
    ylabel('$\mathrm{log}_{10}$ error')
    ha = gca;
    ha.YScale = 'log';

    % legend
    hL_opts = {'Interpreter','latex','FontSize',fontsize_-4,...
        'Location','best','EdgeColor',bcolor};
    legend(hL_opts{:});

end