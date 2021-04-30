%--------------------------------------------------------------------------
% DTQP_MESH_ss_betts.m
% Mesh refinement method for single step methods (based on Betts 4.7)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [T,U,Y,P,F,in,opts] = DTQP_MESH_ss_betts(setup,opts)

% extract
tolerance = opts.dt.meshr.tolerance;
kappa = opts.dt.meshr.errorsafety;
n1 = opts.dt.meshr.ntmaxinterval;
nint = opts.dt.meshr.ntintegration;
maxiters = opts.dt.meshr.maxiters;
ntmax = opts.dt.meshr.ntmax;
storeflag = opts.dt.meshr.storemesh;

% method specific
p = 2; % <- method specific

% extract
displevel = opts.general.displevel;

% change displevel
opts.general.displevel = min(opts.general.displevel,1);

% initialize
doubleflag = false;

% go through each mesh iteration
for iter = 1:maxiters

    % start timer
    tmesh = tic;

    % solve the LQDO problem
    [T,U,Y,P,F,in,opts] = DTQP_multiphase(setup,opts);

    % extract
    nt = in.nt; ny = in.ny;
    A = setup.A; B = setup.B;

    % evaluate time-varying matrices
    At = DTQP_tmultiprod(A,[],T);
    Bt = DTQP_tmultiprod(B,[],T);

    % calculate state derivatives on T
    DY = multiprod(At,Y,[2 3],[2 3]) + multiprod(Bt,U,[2 3],[2 3]);

    %----------------------------------------------------------------------
    % construct continuous representation
    %----------------------------------------------------------------------
    % pchipd for states
    ppYs = pchipd(T',Y',DY');

    % state derivatives using spline derivative
    ppDYs = fnder(ppYs);

    % controls using linear interpolation
    ppUs = griddedInterpolant(T,U,'spline');

    %----------------------------------------------------------------------
    % estimate discretization error
    %----------------------------------------------------------------------
    % absolute local error (4.154)
    eta = zeros(nt-1,ny);
    for k = 1:nt-1
        T2 = linspace(T(k),T(k+1),nint);
        eta(k,:) = trapz(T2,calc_error(T2,ppDYs,ppYs,ppUs,A,B));
    end

    % scale weights (4.157)
    if iter == 1 % only once
        W = max(abs([Y;DY]),[],1);
    end

    % relative local error (4.156)
    Ep = eta./(W+1); % state errors
    ep = max(Ep,[],2); % segment errors

    % maximum error
    epmax = max(ep);

    % (potentially) store results
    if storeflag
        Estorage{iter} = ep; % errors
        Tstorage{iter} = T; % mesh
    end

    % disp iteration
    display_func(displevel,F,max(ep),iter,nt,toc(tmesh),doubleflag)

    % check if mesh error tolerance met
    if epmax < tolerance
        break
    end

    %  average error (4.169)
    epavg = sum(ep)/nt;

    %----------------------------------------------------------------------
    % select primary order for new mesh
    %----------------------------------------------------------------------
    % if the error is equidistributed for the low-order method
    if (p < 4) && (epmax <= 2*epavg)
        p = 4;
    elseif (p < 4) && iter > 2
        p = 4;
    end

    %----------------------------------------------------------------------
    % estimate order reduction
    %----------------------------------------------------------------------
    % compare the current and old grids to compute rk from
    if iter == 1
        r = ones(nt-1,1);
    else

        % replicate indices based on added subintervals
        Ia = repelem(1:length(I),I+1);

        % expand to match eta
        theta2 = theta(Ia,:);
        I2 = I(Ia);

        % (4.160)
        rhat = p + 1 - log(theta2./eta)./log(1+I2);

        % estimated order reduction (4.161)
        rt = max(0,min(rhat,p));

        % determine single value for each segment
        r = round(max(rt,[],2));

    end

    % assign old errors
    theta = eta;

    % create the new mesh
    [opts,I,doubleflag] = create_new_mesh(Ep,ep,nt,tolerance,kappa,n1,p,r,T,opts);

    % check if maximum time points exceeded
    if ntmax < length(opts.dt.t)
        break
    end

end

% return displevel
opts.general.displevel = displevel;

% store mesh
if storeflag
    in.meshr.errors = Estorage;
    in.meshr.T = Tstorage;
    in.meshr.tol = tolerance;
end

end

% calculate errors in the state approximations
function E = calc_error(T,ppDYs,ppYs,ppUs,A,B)

% ensure row vector
T = T(:)';

% calculate state derivatives using spline derivative
DYs = ppval(ppDYs,T);

% calculate states using spline
Y = ppval(ppYs,T)';

% calculate controls using some interpolating method
U = ppUs(T');

% evaluate time-varying matrices
At = DTQP_tmultiprod(A,[],T);
Bt = DTQP_tmultiprod(B,[],T);

% calculate state derivatives on T
DY = multiprod(At,Y,[2 3],[2 3]) + multiprod(Bt,U,[2 3],[2 3]);

% calculate error
E = abs(DYs'-DY);

end

% create the new mesh
function [opts,I,doubleflag] = create_new_mesh(Ep,ep,nt,tol,kappa,n1,p,r,T,opts)

% maximum number of points in the new mesh (old is nt points)
Nmax = 2*nt-1;

% copy original error values
Eporig = Ep;

% initialize array with no points added to each interval
I = zeros(nt-1,1);

% determine where to add new points
while nt < Nmax

    % find interval with maximum error
    [emax,imax] = max(ep);

    % errors is within tolerance
    if emax <= kappa*tol
       break
    end

    % extract current number of points in this interval
    It = I(imax);

    % determine if already too many points in the interval
    if It + 1 > n1
        ep(imax,:) = 0;
        continue
    else % we should add the point
        It = It + 1;

        % add a point to the interval
        I(imax) = It;

        % increment the total number of points
        nt = nt + 1;
    end

    % update the predicted error for interval (4.164)
    ee = Eporig(imax,:)*(1/(1+It)).^(p-r(imax)+1);
    Ep(imax,:) = ee;
    ep(imax) = max(ee);

end

% construct new mesh
Tnew = cell(length(I),1);
for k = 1:length(I)
    Tnew{k} = linspace(T(k),T(k+1),2+I(k));
end

% create the new mesh
Tnew = unique(horzcat(Tnew{:}));

% add the new mesh to opts
opts.dt.mesh = 'USER';
opts.dt.t = Tnew;

% check if the mesh was doubled
doubleflag = nt == Nmax;

end

% local display function for command window
function display_func(displevel,F,eRel,iter,nt,t,doubleflag)

% display to the command window
if (displevel > 1) % verbose

    % handle edge cases
    if isinf(eRel)
        eRelstr = blanks(length(sprintf('%1.2e',1)));
    else
        eRelstr = sprintf('%1.2e',eRel);
    end
    if isinf(F)
        Fstr = blanks(length(sprintf('%1.3e',1)));
    else
        Fstr = sprintf('%1.3e',F);
    end
    if doubleflag
        ntstr = sprintf('*%5i',nt);
    else
        ntstr = sprintf(' %5i',nt);
    end
    tstr = sprintf('%1.2e',t);

    % initial headers
    if iter == 1
        flag = 'line'; DTQP_commandWindowTasks %#ok<NASGU>
        disp('| iter |     nt |         F |    error |    t (s) |')
    end

    % values
    disp(['| ',sprintf(' %3i',iter),...
        ' | ',ntstr,...
        ' | ', Fstr,...
        ' | ', eRelstr,...
        ' | ', tstr,...
        ' | ']);
end
end