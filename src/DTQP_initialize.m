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
if isfield(setup,'p')
    in.p = setup.p;
%     setup = rmfield(setup,'p'); % kind of slow, removing for now
else
    in.p = [];
end

%% create mesh
if ~isfield(setup,'t0') % initial time
    setup.t0 = 0; in.t0 = 0; % default is 0
else
    in.t0 = setup.t0;
end
in.tf = setup.tf; % final time
[t,w,D] = DTQP_pts(in,dt); % mesh, quadrature wights, differentiation matrix
in.nt = length(t); % number of time points
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

% counts
in.ny = max([size(setup.A,1),size(setup.B,1),size(setup.G,1),size(setup.d,1)]); % states
in.nu = size(setup.B,2); % controls
in.np = size(setup.G,2); % parameters
in.nd = size(setup.d,2); % disturbances
in.nx = (in.nu+in.ny)*in.nt + in.np; % optimization variables

% create a zero A matrix when not present
if isempty(setup.A)
	setup.A = zeros(in.ny);
end

%% indices
i = cell(7,1);
i{1} = 1:in.nu; % control indices
i{2} = in.nu+1:(in.nu+in.ny); % state indices
i{3} = (in.nu+in.ny)+1:(in.nu+in.ny)+in.np; % parameter indices
i{4} = i{2}; % initial state indices
i{5} = i{2}; % final state indices
i{6} = i{1}; % initial control indices
i{7} = i{1}; % final control indices
in.i = i;

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

%% add empty fields if not present
if ~isfield(setup,'Y'), setup.Y = []; end
if ~isfield(setup,'Z'), setup.Z = []; end
if ~isfield(setup,'LB'), setup.LB = []; end
if ~isfield(setup,'UB'), setup.UB = []; end
if ~isfield(setup,'LY'), setup.LY = []; end
if ~isfield(setup,'LZ'), setup.LZ = []; end
if ~isfield(setup,'scaling'), setup.scaling = []; end