%--------------------------------------------------------------------------
% DTQP_initialize.m
% Initialize various structures and variables
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [setup,p] = DTQP_initialize(setup,dt)

%% create mesh
p = setup.p;
p.nt = dt.nt;
if ~isfield(p,'t0')
    p.t0 = 0; % default is 0
end
p = DTQP_createT(p,dt);
p.h = diff(p.t);

% calculate midpoints in time grid
if strcmpi(dt.quadrature,'CQHS') || any(strcmpi(dt.defects,{'HS','RK4'}))
    t1 = [0;p.t]; t2 = [p.t;0];
    t3 = t1 + (t2-t1)/2; % midpoints endpoints incorrect
    p.tm = t3(2:end-1);
else
    p.tm = [];
end

%% pseudospectral methods
% differentiation matrix
if strcmpi(dt.defects,'PS')
    tau = ( 2/(p.tf-p.t0)*p.t - (p.tf+p.t0)/(p.tf-p.t0) );
    switch lower(dt.mesh)
        case 'lgl'
            p.D = DTQP_Dmatrix_LGL(tau); % sparse differentiation matrix
        case 'cgl'
            p.D = DTQP_Dmatrix_CGL(tau); % sparse differentiation matrix
    end
else
    p.D = [];
end
% Gaussian quadrature weights
if strcmpi(dt.quadrature,'G') && strcmpi(dt.mesh,'LGL')
	tau = ( 2/(p.tf-p.t0)*p.t - (p.tf+p.t0)/(p.tf-p.t0) );		
    p.w = DTQP_weights_LGL(tau); % sparse Gaussian quadrature weights
elseif strcmpi(dt.quadrature,'CC') && strcmpi(dt.mesh,'CGL')
	tau = ( 2/(p.tf-p.t0)*p.t - (p.tf+p.t0)/(p.tf-p.t0) );		
    p.w = DTQP_weights_CGL(tau); % sparse Gaussian quadrature weights
else
    p.w = [];
end

%% counts
if ~isfield(setup,'A'), setup.A = []; end % needs to be zero
if ~isfield(setup,'B'), setup.B = []; end
if ~isfield(setup,'G'), setup.G = []; end
if ~isfield(setup,'d'), setup.d = []; end

if ~isfield(setup,'L'), setup.L = []; end
if ~isfield(setup,'M'), setup.M = []; end
if ~isfield(setup,'c'), setup.c = []; end

p.ns = max([size(setup.A,1),size(setup.B,1),size(setup.G,1),size(setup.d,1)]); % number of states
p.nu = size(setup.B,2); % number of controls
p.np = size(setup.G,2); % number of parameters
p.nd = size(setup.d,2); % number of disturbances

p.nx = (p.nu+p.ns)*p.nt + p.np; % number of optimization variables

if isempty(setup.A)
	setup.A = zeros(p.ns);
end

%% add fields if not present
% go through Lagrange terms
Ltemp = setup.L;
setup = rmfield(setup,'L');
for i = 1:length(Ltemp)
    % check if there is a nonzero left field
    Lflag = 0;
    if isfield(Ltemp(i),'left')
        if Ltemp(i).left ~= 0
            Lflag = 1;
        end
    end

    % check if there is a nonzero right field
    Rflag = 0;
    if isfield(Ltemp(i),'right')
        if Ltemp(i).right ~= 0
            Rflag = 1;
        end
    end
    
    % determine the order
    if (Lflag == 0) && (Rflag == 0) % constant term
        if isfield(setup,'cL')
            N = length(setup.cL);
        else
            N = 0;
        end
        setup.cL(N+1).matrix = Ltemp(i).matrix;
    elseif (Lflag == 0) && (Rflag == 1) % linear term
        if isfield(setup,'l')
            N = length(setup.l);
        else
            N = 0;
        end
        setup.l(N+1).right = Ltemp(i).right;
        setup.l(N+1).matrix = Ltemp(i).matrix;
    elseif (Lflag == 1) && (Rflag == 0) % linear term (improper ordering)
        if isfield(setup,'l')
            N = length(setup.l);
        else
            N = 0;
        end
        setup.l(N+1).right = Ltemp(i).left; % move left to right
        setup.l(N+1).matrix = Ltemp(i).matrix'; % transpose
    else % quadratic term
        if isfield(setup,'L')
            N = length(setup.L);
        else
            N = 0;
        end
        setup.L(N+1).left = Ltemp(i).left;
        setup.L(N+1).right = Ltemp(i).right;
        setup.L(N+1).matrix = Ltemp(i).matrix;
    end
    
end

% go through Mayer terms
Mtemp = setup.M;
setup = rmfield(setup,'M');
for i = 1:length(Mtemp)
    % check if there is a nonzero left field
    Lflag = 0;
    if isfield(Mtemp(i),'left')
        if Mtemp(i).left ~= 0
            Lflag = 1;
        end
    end

    % check if there is a nonzero right field
    Rflag = 0;
    if isfield(Mtemp(i),'right')
        if Mtemp(i).right ~= 0
            Rflag = 1;
        end
    end

    % determine the order
    if (Lflag == 0) && (Rflag == 0) % constant term
        if isfield(setup.c,'M')
            N = length(setup.cM);
        else
            N = 0;
        end
        setup.cM(N+1) = Mtemp(i).matrix;
    elseif (Lflag == 0) && (Rflag == 1) % linear term
        if isfield(setup,'m')
            N = length(setup.m);
        else
            N = 0;
        end
        setup.m(N+1).right = Mtemp(i).right;
        setup.m(N+1).matrix = Mtemp(i).matrix;
    else % quadratic term
        if isfield(setup,'M')
            N = length(setup.m);
        else
            N = 0;
        end
        setup.M(N+1).left = Mtemp(i).left;
        setup.M(N+1).right = Mtemp(i).right;
        setup.M(N+1).matrix = Mtemp(i).matrix;
    end
    
end

% add empty fields if not present
if ~isfield(setup,'L'), setup.L = []; end
if ~isfield(setup,'M'), setup.M = []; end
if ~isfield(setup,'c'), setup.c = []; end
if ~isfield(setup,'l'), setup.l = []; end
if ~isfield(setup,'m'), setup.m = []; end
if ~isfield(setup,'cL'), setup.cL = []; end
if ~isfield(setup,'cM'), setup.cM = []; end
if ~isfield(setup,'Y'), setup.Y = []; end
if ~isfield(setup,'Z'), setup.Z = []; end
if ~isfield(setup,'LB'), setup.LB = []; end
if ~isfield(setup,'UB'), setup.UB = []; end
if ~isfield(setup,'scaling'), setup.scaling = []; end
if ~isfield(setup,'LY'), setup.LY = []; end
if ~isfield(setup,'LZ'), setup.LZ = []; end

%% indices
p.i{1} = 1:p.nu; % control indices
p.i{2} = p.nu+1:(p.nu+p.ns); % state indices
p.i{3} = (p.nu+p.ns)+1:(p.nu+p.ns)+p.np; % parameter indices
p.i{4} = p.i{2}; % initial state indices
p.i{5} = p.i{2}; % final state indices
p.i{6} = p.i{1}; % initial control indices
p.i{7} = p.i{1}; % final control indices