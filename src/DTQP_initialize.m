%--------------------------------------------------------------------------
% DTDP_initialize.m
% Initialize various structures and variables
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber), University of 
% Illinois at Urbana-Champaign
% Project link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------

%% create mesh
p = setup.p;
p = DTQP_createT(p,opts);
p.h = diff(p.t);

% calculate midpoints in time grid
if strcmp(opts.Quadmethod,'CQHS') || any(strcmp(opts.Defectmethod,{'HS','RK4'}))
    t1 = [0;p.t]; t2 = [p.t;0];
    t3 = t1 + (t2-t1)/2; % midpoints endpoints incorrect
    p.tm = t3(2:end-1);
end

%% pseudospectral methods
% differentiation matrix
if strcmp(opts.Defectmethod,'PS') || strcmp(opts.Defectmethod,'PS_old')
    tau = ( 2/(p.tf-p.t0)*p.t - (p.tf+p.t0)/(p.tf-p.t0) );
    switch opts.NType
        case 'LGL'
            p.D = DTQP_Dmatrix_LGL(tau); % sparse differentiation matrix
        case 'CGL'
            p.D = DTQP_Dmatrix_CGL(tau); % sparse differentiation matrix
    end
end
% Gaussian quadrature weights
if strcmp(opts.Quadmethod,'G') && strcmp(opts.NType,'LGL')
	tau = ( 2/(p.tf-p.t0)*p.t - (p.tf+p.t0)/(p.tf-p.t0) );		
    p.w = DTQP_weights_LGL(tau); % sparse Gaussian quadrature weights
end
if strcmp(opts.Quadmethod,'CC') && strcmp(opts.NType,'CGL')
	tau = ( 2/(p.tf-p.t0)*p.t - (p.tf+p.t0)/(p.tf-p.t0) );		
    p.w = DTQP_weights_CGL(tau); % sparse Gaussian quadrature weights
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

p.nt = length(p.t);
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
        if isfield(setup.c,'L')
            N = length(setup.c.L);
        else
            N = 0;
        end
        setup.c.L(N+1) = Ltemp(i).matrix;
    elseif (Lflag == 0) && (Rflag == 1) % linear term
        if isfield(setup,'l')
            N = length(setup.l);
        else
            N = 0;
        end
        setup.l(N+1).right = Ltemp(i).right;
        setup.l(N+1).matrix = Ltemp(i).matrix;
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
            N = length(setup.c.M);
        else
            N = 0;
        end
        setup.c.M(N+1) = Mtemp(i).matrix;
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
if ~isfield(setup.c,'L'), setup.c.L = []; end
if ~isfield(setup.c,'M'), setup.c.M = []; end
if ~isfield(setup,'Y'), setup.Y = []; end
if ~isfield(setup,'Z'), setup.Z = []; end
if ~isfield(setup,'LB'), setup.LB = []; end
if ~isfield(setup,'UB'), setup.UB = []; end
if ~isfield(setup,'scaling'), setup.scaling = []; end

%% indices
p.i{1} = 1:p.nu; % control indices
p.i{2} = p.nu+1:(p.nu+p.ns); % state indices
p.i{3} = (p.nu+p.ns)+1:(p.nu+p.ns)+p.np; % parameter indices
p.i{4} = p.i{2}; % initial state indices
p.i{5} = p.i{2}; % final state indices