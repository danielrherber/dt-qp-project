%--------------------------------------------------------------------------
% DTQP_qlin_plots.m
% Create the plots for the quasilinearization iterations
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function DTQP_qlin_plots(T,Y,U,P,iter)

% figure numbers
hfY = 1; hfU = 2; hfP = 3;

% initial actions
if (iter==1)
    close all
    if ~isempty(Y)
        hf = figure(hfY); hold on; hf.Color = 'w'; title('States'); xlabel("$t$");
    end
    if ~isempty(U)
        hf = figure(hfU); hold on; hf.Color = 'w'; title('Controls'); xlabel("$t$");
    end
    if ~isempty(P)
        hf = figure(hfP); hold on; hf.Color = 'w'; title('Parameters'); xlabel("iteration");
    end
end

% determine current line color
switch iter
    case 0 % final iteration?
        c = 'r';
    case 1 % initial iteration?
        c = 'k';
    otherwise % intermediate iteration
        c = [0.5 0.5 0.5];
end

% plot states
if ~isempty(Y)
    figure(hfY)
    plot(T,Y,'color',c,'linewidth',2);
end

% plot controls
if ~isempty(U)
    figure(hfU)
    plot(T,U,'color',c,'linewidth',2);
end

% plot parameters
if ~isempty(P)
    figure(hfP)
    plot(iter,P,'.','color',c);
end

end