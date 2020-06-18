%------------------------------------------------------------------------------
% DTQP_default_opts.m
% Obtain default options
%------------------------------------------------------------------------------
%
%------------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%------------------------------------------------------------------------------
function [setup,opts] = DTQP_default_opts(setup,opts)

%--------------------------------------------------------------------------
% START: general options
%--------------------------------------------------------------------------
% initialize general options structure
if ~isfield(opts,'general')
    opts.general = [];
end

% plots
if ~isfield(opts.general,'plotflag')
    opts.general.plotflag = 1; % create the plots
    % opts.general.plotflag = 0; % don't create the plots
end

% save the solution and plots to disk (custom plotting function)
if ~isfield(opts.general,'saveflag')
    % opts.general.saveflag = 1; % save
    opts.general.saveflag = 0; % don't save
end

% name of the example
if ~isfield(opts.general,'name')
    % opts.general.name = mfilename;
    opts.general.name = 'DTQP_Example';
end

% path for saving
if ~isfield(opts.general,'path')
    opts.general.path = mfoldername(mfilename('fullpath'),'_private');
end

% controls the displayed diagnostics in the command window
if ~isfield(opts.general,'displevel')
    opts.general.displevel = 2; % verbose
    % opts.general.displevel = 1; % minimal
    % opts.general.displevel = 0; % none
end

% scaling to the constraint rows
if ~isfield(opts.general,'scalerowflag')
    opts.general.scalerowflag = true; % enabled
    % opts.general.scalerowflag = false; % disabled
end
%--------------------------------------------------------------------------
% END: general options
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% START: direct transcription specific
%--------------------------------------------------------------------------
% initialize direct transcription-specific options structure
if ~isfield(opts,'dt')
    opts.dt(1).refinement = []; % initialize
end

%--------------------------------------------------------------------------
% mesh refinement algorithm (first phase)
defaultflag = 0; % initialize
if ~isfield(opts.dt(1),'meshr') % not present
    opts.dt(1).meshr = DTQP_meshr_default_opts([]); % set as the default
    defaultflag = 1;
elseif isempty(opts.dt(1).meshr) % empty
    opts.dt(1).meshr = DTQP_meshr_default_opts([]); % set as the default
    defaultflag = 1;
else
    opts.dt(1).meshr = DTQP_meshr_default_opts(opts.dt(1).meshr); %
end
if defaultflag && (opts.general.displevel > 1) % minimal
    disp(['using default mesh refinement ',opts.dt(1).meshr.method])
end

%--------------------------------------------------------------------------
% defect constraint method (first phase)

% default = 'ZO'; % zero-order hold
% default = 'EF'; % Euler forward
% default = 'Huen'; % Huen's method
% default = 'ModEF'; % modified Euler method
default = 'TR'; % trapezoidal
% default = 'HS'; % Hermite-Simpson
% default = 'RK4'; % fourth-order Runge-Kutta
% default = 'PS'; % pseudospectral (both LGL and CGL meshes)

defaultflag = 0; % initialize
if ~isfield(opts.dt(1),'defects') % not present
    opts.dt(1).defects = default; % set as the default
    defaultflag = 1;
elseif isempty(opts.dt(1).defects) % empty
    opts.dt(1).defects = default; % set as the default
    defaultflag = 1;
end
if defaultflag &&(opts.general.displevel > 1) % minimal
    disp(['using default defect constraint method ',opts.dt(1).defects])
end

%--------------------------------------------------------------------------
% quadrature method (first phase)

% default = 'CEF'; % composite Euler forward
default = 'CTR'; % composite trapezoidal
% default = 'CQHS'; % composite quadratic Hermite-Simpson
% default = 'G'; % Gaussian
% default = 'CC'; % Clenshaw-Curtis

defaultflag = 0; % initialize
if ~isfield(opts.dt(1),'quadrature') % not present
    opts.dt(1).quadrature = default; % set as the default
    defaultflag = 1;
elseif isempty(opts.dt(1).quadrature) % empty
    opts.dt(1).quadrature = default; % set as the default
    defaultflag = 1;
end
if defaultflag &&(opts.general.displevel > 1) % minimal
    disp(['using default quadrature method ',opts.dt(1).quadrature])
end

%--------------------------------------------------------------------------
% mesh type (first phase)

default = 'ED'; % equidistant nodes
% default = 'LGL'; % Legendre-Gauss-Lobatto nodes
% default = 'CGL'; % Chebyshev-Gauss-Lobatto nodes
% default = 'USER'; % user-defined nodes

defaultflag = 0; % initialize
if ~isfield(opts.dt(1),'mesh') % not present
    opts.dt(1).mesh = default; % set as the default
    defaultflag = 1;
elseif isempty(opts.dt(1).mesh) % empty
    opts.dt(1).mesh = default; % set as the default
    defaultflag = 1;
end
if defaultflag &&(opts.general.displevel > 1) % minimal
        disp(['using default mesh type ',opts.dt(1).mesh])
end

%--------------------------------------------------------------------------
% number of nodes (first phase)

default = 100; % 100 nodes

if isfield(opts.dt(1).meshr,'method')
    if ~strcmpi(opts.dt(1).meshr.method,'none')
        default = nan; % will be set in DTQP_meshr_default_opts.m
    end
end

defaultflag = 0; % initialize
if ~isfield(opts.dt(1),'nt') % not present
    opts.dt(1).nt = default; % set as the default
    defaultflag = 1;
elseif isempty(opts.dt(1).nt) % empty
    opts.dt(1).nt = default; % set as the default
    defaultflag = 1;
end
if defaultflag &&(opts.general.displevel > 1) % minimal
        disp(['using default number of nodes ',opts.dt(1).nt])
end

%--------------------------------------------------------------------------
% phase specific
nphase = length(setup); % number of phases
for phs = 2:nphase

    % potentially copy phase information
    if phs > length(opts.dt)
        DT = [];
    else
        DT = opts.dt(phs);
    end

    % mesh refinement algorithm (same for all phases)
    if ~isfield(DT,'meshr')
        DT.meshr = opts.dt(1).meshr; % from first phase
    elseif isempty(DT.meshr)
        DT.meshr = opts.dt(1).meshr; % from first phase
    else
        DT.meshr = DTQP_meshr_default_opts(DT.meshr); %
    end

    % defect constraint method
    if ~isfield(DT,'defects')
        DT.defects = opts.dt(1).defects; % from first phase
    elseif isempty(DT.defects)
        DT.defects = opts.dt(1).defects; % from first phase
    end

    % quadrature method
    if ~isfield(DT,'quadrature')
        DT.quadrature = opts.dt(1).quadrature; % from first phase
    elseif isempty(DT.quadrature)
        DT.quadrature = opts.dt(1).quadrature; % from first phase
    end

    % mesh type
    if ~isfield(DT,'mesh')
        DT.mesh = opts.dt(1).mesh; % from first phase
    elseif isempty(DT.mesh)
        DT.mesh = opts.dt(1).mesh; % from first phase
    end

    % number of nodes
    if ~isfield(DT,'nt')
        DT.nt = opts.dt(1).nt; % from first phase
    elseif isempty(DT.nt)
        DT.nt = opts.dt(1).nt; % from first phase
    end

    % assign and clear
    opts.dt(phs) = DT;
    clear DT

end
%--------------------------------------------------------------------------
% END: direct transcription specific
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% START: quadratic programming specific
%--------------------------------------------------------------------------
% initialize quadratic programming-specific options structure
if ~isfield(opts,'qp')
    opts.qp = [];
end

% reordering of the optimization variables (see DTQP_reorder.m)
if ~isfield(opts.qp,'reorder')
    opts.qp.reorder = 0; % don't reorder
    % opts.qp.reorder = 1; % reorder
end

% solver
if ~isfield(opts.qp,'solver')
    opts.qp.solver = 'quadprog'; % MATLAB quadprog
    % opts.qp.solver = 'cvx'; % see DTQP_solver_cvx.m
    % opts.qp.solver = 'qpoases'; % see DTQP_solver_qpoases.m
    % opts.qp.solver = 'qpip'; % (to be added)
    % opts.qp.solver = 'ooqp'; % (to be added)
end

% get default options for the selected solver
opts = DTQP_solver_default_opts(opts);
%--------------------------------------------------------------------------
% END: quadratic programming specific
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% START: other
%--------------------------------------------------------------------------

% get default options for NLDO problem
opts = DTQP_nonlin_defaultOpts(opts);

%--------------------------------------------------------------------------
% END: other
%--------------------------------------------------------------------------
end