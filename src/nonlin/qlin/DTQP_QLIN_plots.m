%--------------------------------------------------------------------------
% DTQP_QLIN_plots.m
% Create the plots for the quasilinearization iterations
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP_QLIN_plots(T,Y,U,P,iter)

% determine what design variable types are present
n = 0;
if ~isempty(Y)
    n = n + 1;
    IY = n;
end
if ~isempty(U)
    n = n + 1;
    IU = n;
end
if ~isempty(P)
    n = n + 1;
    IP = n;
end

% initial actions
if (iter==1)
    flag = 'preliminary'; DTQP_plotCommon; %#ok<NASGU>
    hf = figure(1); hold on; hf.Color = wcolor;
    hf.Position(3:4) = [1200 420]; hf.Position(1) = 500;
    if ~isempty(Y)
        subplot(1,n,IY,'align'); hold on; xlabel("$t$");
        title('States','fontsize',fontsize_+2);
        flag = 'axis'; DTQP_plotCommon; %#ok<NASGU>
        legend('location','best');
    end
    if ~isempty(U)
        subplot(1,n,IU,'align'); hold on; hf.Color = 'w'; xlabel("$t$");
        title('Controls','fontsize',fontsize_+2);
        flag = 'axis'; DTQP_plotCommon; %#ok<NASGU>
        legend('location','best');
    end
    if ~isempty(P)
        subplot(1,n,IP,'align'); hold on; hf.Color = 'w'; xlabel("iteration");
        title('Parameters','fontsize',fontsize_+2);
        flag = 'axis'; DTQP_plotCommon; %#ok<NASGU>
        legend('location','best');
    end
end

% determine current line color
if iter == 1 % initial iteration?
    c = 'k';
elseif iter < 0 % final iteration?
    c = 'r';
else % intermediate iteration
    c = [0.5 0.5 0.5];
end

% plot states
if ~isempty(Y)
    ha = subplot(1,n,IY);
    plot(T,Y,'color',c,'linewidth',2);
    modLegend(ha,iter);
end

% plot controls
if ~isempty(U)
    ha = subplot(1,n,IU);
    plot(T,U,'color',c,'linewidth',2);
    modLegend(ha,iter);
end

% plot parameters
if ~isempty(P)
    ha = subplot(1,n,IP);
    plot(abs(iter),P,'.','color',c);
    modLegend(ha,iter);
end

% draw the current iteration
drawnow;

end

function modLegend(ha,iter)

if iter == 1
    ha.Legend.String{1} = 'Initial';
    ha.Legend.String(2:end) = [];
elseif iter < 0
    ha.Legend.String{3} = 'Final';
    ha.Legend.String(4:end) = [];
elseif iter == 2
    ha.Legend.String{2} = 'Intermediate';
    ha.Legend.String(3:end) = [];
else
    ha.Legend.String(3:end) = [];
end

end