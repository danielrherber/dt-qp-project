%--------------------------------------------------------------------------
% DTQP_initialize.m
% Initialize various structures and variables
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup,in] = DTQP_initialize(setup,dt)

% set user parameter structure
if isfield(setup,'auxdata')
    in.auxdata = setup.auxdata;
%     setup = rmfield(setup,'auxdata'); % kind of slow, removing for now
else
    in.auxdata = [];
end

%% create mesh
if ~isfield(setup,'t0') % initial time
    setup.t0 = 0; in.t0 = 0; % default is 0
else
    in.t0 = setup.t0;
end
in.tf = setup.tf; % final time
[t,w,D] = DTQP_MESH_pts(in,dt); % mesh, quadrature wights, differentiation matrix
nt = length(t); in.nt = nt; % number of time points
in.t = t; in.w = w(:); in.D = D;
h = diff(t); in.h = h; % time steps

% setup = rmfield(setup,'t0'); % kind of slow, removing for now
% setup = rmfield(setup,'tf'); % kind of slow, removing for now

% calculate required interior points in time mesh
if strcmpi(dt.quadrature,'CQHS') || any(strcmpi(dt.defects,{'HS','RK4'}))
    in.tm = t(1:end-1) + h/2; % midpoints
elseif strcmpi(dt.defects,'Huen')
    in.tm = t(1:end-1) + (2/3)*h;
else
    in.tm = [];
end

%% counts
% initialize dynamic matrices if not present
if ~isfield(setup,'A'), setup.A = []; end
if ~isfield(setup,'B'), setup.B = []; end
if ~isfield(setup,'G'), setup.G = []; end
if ~isfield(setup,'d'), setup.d = []; end

% field containing number of controls, states, and parameters
if ~isfield(setup,'n')
    ns = [];
else
    ns = setup.n;
end

% states
if isfield(ns,'ny') && ~isempty(ns.ny)
    in.ny = ns.ny;
else
    in.ny = max([size(setup.A,1),size(setup.B,1),size(setup.G,1),size(setup.d,1)]);
end
ny = in.ny;

% controls
if isfield(ns,'nu') && ~isempty(ns.nu)
    in.nu = ns.nu;
else
    in.nu = size(setup.B,2);
end
nu = in.nu;

% parameters
if isfield(ns,'np') && ~isempty(ns.np)
    in.np = ns.np;
else
    in.np = size(setup.G,2);
end
np = in.np;

% disturbances
in.nd = size(setup.d,2);

% optimization variables
in.nx = (nu+ny)*nt + np;

% create a zero A matrix when not present
if isempty(setup.A)
	setup.A = zeros(ny);
end

%% indices
% infinite-dimensional problem
i = cell(7,1);
i{1} = 1:nu; % control indices
i{2} = nu+1:(nu+ny); % state indices
i{3} = (nu+ny)+1:(nu+ny)+np; % parameter indices
i{4} = i{2}; % initial state indices
i{5} = i{2}; % final state indices
i{6} = i{1}; % initial control indices
i{7} = i{1}; % final control indices
in.i = i;

% finite-dimensional problem
in.I_stored = [reshape(1:(nu+ny)*nt,nt,[]),repmat((nu+ny)*nt+(1:np),nt,1)];

%% go through Lagrange terms
% add fields if not present
if isfield(setup,'L')
    Ltemp = setup.L;
else
    Ltemp = [];
end

% initialize structures
Lfields = {'left','right','matrix';{},{},{}};
L = struct(Lfields{:});
lfields = {'right','matrix';{},{}};
l = struct(lfields{:});
cLfields = {'matrix';{}};
cL = struct(cLfields{:});

% go through each entry in L
for k = 1:length(Ltemp)

    % check if there is a nonzero left field
    Lflag = 0;
    if isfield(Ltemp(k),'left')
        if Ltemp(k).left ~= 0
            Lflag = 1;
        end
    end

    % check if there is a nonzero right field
    Rflag = 0;
    if isfield(Ltemp(k),'right')
        if Ltemp(k).right ~= 0
            Rflag = 1;
        end
    end

    % determine correct structure
    if (Lflag == 0) && (Rflag == 0) % constant term
        cL(end+1).left = 0;
        cL(end).right = 0;
        cL(end).matrix = Ltemp(k).matrix;
    elseif (Lflag == 0) && (Rflag == 1) % linear term
        l(end+1).left = 0;
        l(end).right = Ltemp(k).right;
        l(end).matrix = Ltemp(k).matrix;
    elseif (Lflag == 1) && (Rflag == 0) % linear term (improper ordering)
        l(end+1).right = Ltemp(k).left; % move left to right
        l(end).left = 0;
        l(end).matrix = Ltemp(k).matrix'; % transpose
    else % quadratic term
        L(end+1) = Ltemp(k); % copy all fields
    end

end

% assign to setup structure
setup.L = L; setup.l = l; setup.cL = cL;

%% go through Mayer terms
% add fields if not present
if isfield(setup,'M')
    Mtemp = setup.M;
else
    Mtemp = [];
end

% initialize structures
Mfields = {'left','right','matrix';{},{},{}};
M = struct(Mfields{:});
mfields = {'right','matrix';{},{}};
m = struct(mfields{:});
cMfields = {'matrix';{}};
cM = struct(cMfields{:});

% go through each entry in M
for k = 1:length(Mtemp)

    % check if there is a nonzero left field
    Lflag = 0;
    if isfield(Mtemp(k),'left')
        if Mtemp(k).left ~= 0
            Lflag = 1;
        end
    end

    % check if there is a nonzero right field
    Rflag = 0;
    if isfield(Mtemp(k),'right')
        if Mtemp(k).right ~= 0
            Rflag = 1;
        end
    end

    % determine correct structure
    if (Lflag == 0) && (Rflag == 0) % constant term
        cM(end+1).left = 0;
        cM(end).right = 0;
        cM(end).matrix = Mtemp(k).matrix;
    elseif (Lflag == 0) && (Rflag == 1) % linear term
        m(end+1).left = 0;
        m(end).right = Mtemp(k).right;
        m(end).matrix = Mtemp(k).matrix;
    else % quadratic term
        M(end+1) = Mtemp(k); % copy all fields
    end

end

% assign to setup structure
setup.M = M; setup.m = m; setup.cM = cM;

%% go through scaling terms
% add fields if not present
if ~isfield(setup,'scaling')
    setup.scaling = [];
    scaling = [];
else
    scaling = setup.scaling;
end

% initialize
stemp = struct('right',{},'matrix',{},'constant',{});

% go through each entry in scaling
for k = 1:length(scaling)

    % extract
    s = scaling(k);

    % determine number of variables
    switch s.right
        case 1
            n = nu;
        case 2
            n = ny;
        case 3
            n = np;
    end

    % check matrix field
    if ~isfield(s,'matrix')
        s.matrix = ones(1,n);
    elseif isempty(s.matrix)
        s.matrix = ones(1,n);
    end

    % check constant field
    if ~isfield(s,'constant')
        s.constant = zeros(1,n);
    elseif isempty(s.constant)
        s.constant = zeros(1,n);
    end

    % assign
    stemp(k) = s;

end

% assign to setup structure
setup.scaling = stemp;

%% add empty fields if not present
if ~isfield(setup,'Y'), setup.Y = []; end
if ~isfield(setup,'Z'), setup.Z = []; end
if ~isfield(setup,'LB'), setup.LB = []; end
if ~isfield(setup,'UB'), setup.UB = []; end
if ~isfield(setup,'LY'), setup.LY = []; end
if ~isfield(setup,'LZ'), setup.LZ = []; end