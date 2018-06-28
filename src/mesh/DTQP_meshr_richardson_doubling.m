%--------------------------------------------------------------------------
% DTQP_meshr_richardson_doubling.m
% Simple mesh refinement scheme based on Richardson extrapolation and
% approximately doubling the number of time points
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_meshr_richardson_doubling(setup,opts,solvefun)

% extract options
dt = opts.dt;
if isnan(dt(1).nt)
    ntinit = dt(1).meshr.ntinit; % initial number of steps
else
    ntinit = dt(1).nt;
end
ntmax = dt(1).meshr.ntmax; % max number of time points
etol = dt(1).meshr.tolerance; % relative function tolerance
t = 2; % refinement factor

% initialize
QPcreatetime = 0; QPsolvetime = 0; % timers
eRel = Inf; % initial error
nt = 0; % initial number of time points
iter = 0; % iteration counter
k0 = nan; % convergence rate
fminbndopts = optimset('Display','none','TolX',1e-10); % fminbnd options

% storage elements
ms = []; % scale factors
k0s = []; % k values
Fs = []; % objective function values
Rs = []; % Richardson extrapolation objective function values

% local display flag
displevel = opts.general.displevel;
opts.general.displevel = displevel > 0;

% approximately double the number of time points until convergence
while (eRel > etol) && (nt < ntmax)
    % multiplicative factor
    m = t^iter; ms(end+1) = m; % store
    
    % calculate the new number of time points
    nt = round(m*ntinit);

    % update number of time points
    opts.dt.nt = nt;

    % try to solve the problem with the specified mesh
    try
        % solve
        [T,U,Y,P,F,in,opts] = solvefun(setup,opts);

        % combine timers (NOTE: this currently misses timers with errors)
        if (opts.general.displevel > 0) % minimal
            QPcreatetime = QPcreatetime + in.QPcreatetime;
            QPsolvetime = QPsolvetime + in.QPsolvetime;
        end
    catch
        % reset saved quantities when QP problem fails
        ms = []; Fs = []; k0s = []; F = inf;
    end

    % check if the QP problem did not converge
    if isnan(F) ||isinf(F)
        % reset saved quantities when QP problem did not converge
        ms = []; Fs = []; k0s = []; F = inf;
    else
        % store objective function value
        Fs(end+1) = F;
    end

    if length(Fs) > 2
        k0max = 10;
        % determine the current value of k0 using previous 3 values
        k0 = fminbnd(@(k0) rate_error(k0,Fs(end-2:end),ms(end-2:end)),0,k0max,fminbndopts);

        % leading order errors in A1s
        Ks = k0*(1:length(Fs)); % assumption on the powers

        % Richardson extrapolation for F using previous objective values
        R = richardson(Fs,Ks,t);

        % relative error using Richardson extrapolation
        eRel = abs(1 - F/R); 

        % store values
        k0s(end+1) = k0;
        Rs(end+1) = R;
    elseif length(Fs) == 2
        eRel = abs(1 - Fs(end-1)/F); % relative error with previous iteration
    elseif length(Fs) == 1
        eRel = Inf;
    end

    % display to the command window
    if (displevel > 1) % verbose
        % handle edge cases
        if isnan(k0)
            k0str = blanks(length(sprintf('%1.3e',1)));
        else
            k0str = sprintf('%1.3e',k0);
        end
        if isinf(eRel)
            eRelstr = blanks(length(sprintf('%1.3e',1)));
        else
            eRelstr = sprintf('%1.3e',eRel);
        end
        if isinf(F)
            Fstr = blanks(length(sprintf('%1.3e',1)));
        else
            Fstr = sprintf('%1.3e',F);
        end

        % initial headers
        if iter == 0
            disp('----------------------------------------------------')
            disp('| iter |    nt |         F | est_error |    est_k0 |')
        end

        % values
        disp(['| ',sprintf(' %3i',iter),...
            ' | ', sprintf('%5i',nt),...
            ' | ', Fstr,...
            ' | ', eRelstr,...
            ' | ', k0str,...
            ' | ']); % debugging, remove

        % % required nt
        % ntreq = nt*(eRel/etol)^(1/k0);
    end

	% increment counter
    iter = iter + 1;
end

% end the timer
if (displevel > 0) % minimal
    in.QPtotaltime = toc;
    in.QPcreatetime = QPcreatetime;
    in.QPsolvetime = QPsolvetime;
end

% display to the command window
if (displevel > 1) % verbose
    disp('----------------------------------------------------')
    disp(['QP total time: ', num2str(in.QPtotaltime), ' s'])
end

% reset display level
opts.general.displevel = displevel;

end
% convergence rate error function
function z = rate_error(k0,f,m)
% multiplicative factors
m2 = m(2)/m(1);
m3 = m(3)/m(1);

% error
e = (f(2) + (f(2)-f(1))/(m2^k0-1)) - (f(3) + (f(3)-f(1))/(m3^k0-1));

% square of the error
z = e^2;
end
% Richardson extrapolation
function R = richardson(As,ks,t)
% general recurrence relation for Richardson extrapolation
f = @(Ah,Aht,k,t) (t^k*Aht - Ah)./(t^k - 1);

% compute the triangular extrapolation table
for i = 1:(size(As,2)-1)
    for j = 1:(size(As,2)-i)
        As(:,j) = f(As(:,j),As(:,j+1),ks(i),t);
    end
end

% output extrapolated value
R = As(:,1);
end