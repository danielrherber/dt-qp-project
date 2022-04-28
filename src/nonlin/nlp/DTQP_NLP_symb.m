%--------------------------------------------------------------------------
% DTQP_NLP_symb.m
% Construct MATLAB function and its derivatives for use with ipfmincon
% method
%--------------------------------------------------------------------------
% f must be a string that contains only:
% - y (e.g. y1, y2,...) with ny states
% - u (e.g. u1, u2,...) with nu inputs
% - p (e.g. p1, p2,...) with np parameters
% - additional symbolic variables in varargin{1}
% n is the order of the Taylor polynomial approximation
%--------------------------------------------------------------------------
% Type 1: just a linear time-varying vector of equations
% Type 2: just a quadratic time-varying vector of equations
% Type 3: linearization of a nonlinear vector of equations
% Type 4: quadraticization of a nonlinear scalar equation
%--------------------------------------------------------------------------
% Primary contributor: Daniel R. Herber (danielrherber on GitHub)
% Link: https://github.com/danielrherber/dt-qp-project
%--------------------------------------------------------------------------
function [E,opts] = DTQP_NLP_symb(f,in,linflag,pathboundaryflag,opts)

% (potentially) start the timer
if (opts.general.displevel > 0) % minimal
    opts.timer.t2 = tic; % start timer
end

%--------------------------------------------------------------------------
% parse inputs
%--------------------------------------------------------------------------
% number of states
if isfield(in,'ny')
	ny = in.ny;
else
    ny = 0;
end

% number of controls
if isfield(in,'nu')
	nu = in.nu;
else
    nu = 0;
end

% number of parameters
if isfield(in,'np')
	np = in.np;
else
    np = 0;
end

% number of optimization variables
nx = nu + ny + np;

% extra symbolic parameters
if isfield(in,'parameter_list')
	param = in.parameter_list;
else
    param = [];
end

output = 1; % vectorized matlab function
% output = 2; % interp1-compatible matlab function
% output = 3; % symbolic functions

%--------------------------------------------------------------------------
% initialize the symbolic variables
%--------------------------------------------------------------------------
strY = []; strU = []; strP = []; strYi = []; strYf = [];

% add controls
for idx = 1:nu
    strU = [strU,'u',num2str(idx),' '];
end

% add states
for idx = 1:ny
    strY = [strY,'y',num2str(idx),' '];
end

% add parameters
for idx = 1:np
    strP = [strP,'p',num2str(idx),' '];
end

% add initial states
for idx = 1:ny
    strYi = [strYi,'yi',num2str(idx),' '];
end

% add final states
for idx = 1:ny
    strYf = [strYf,'yf',num2str(idx),' '];
end

% create symbolic variables
eval(['syms t dummy123 ',strU,strY,strP,strYi,strYf,param,' real'])

eval(['U = [',strU,'];'])
eval(['Y = [',strY,'];'])
eval(['P = [',strP,'];'])
eval(['Yi = [',strYi,'];'])
eval(['Yf = [',strYf,'];'])

if isempty(param)
    PARAM = dummy123;
else
	eval(['PARAM = [',param,'];'])
end

%--------------------------------------------------------------------------
% create the vectorized functions
%--------------------------------------------------------------------------
% create symbolic object
F = eval(f);

% ensure column vector
F = reshape(F,[],1);

% vector of optimization variables
X = [U,Y,P,Yi,Yf];

% collection of t, parameters, and optimization variables
in1 = {t,PARAM,X};

% return only function unless symbolic method
if ~strcmpi(opts.method.derivatives,'symbolic')

    % vectorized matlab function
    E.f = sym2matrixfun(F,in1,output);

    % treat all functions as nonlinear
    E.Ilin = [];
    E.Inon = 1:length(F);

    % (potentially) end the timer
    if (opts.general.displevel > 0) % minimal
        opts.timer.sym = opts.timer.sym + toc(opts.timer.t2); % start timer
    end

    return
end

% Jacobian (first partial derivatives)
DF = jacobian(F,X');

% check for the linear equations
if linflag

    % initialize all as nonlinear
    Ilin = false(1,size(DF,1));

    % go through each equation
    for k = 1:size(DF,1)

        % check if the equation is linear with respect to X
        try
            equationsToMatrix(F(k,:),X);
            Ilin(k) = true;
        catch

        end
    end

    % extract linear equations
    Flin = F(Ilin);

    % remove linear equations of F and DF
    F(Ilin) = [];
    DF(Ilin,:) = [];

    % convert linear equations to matrix form
    [L,O] = equationsToMatrix(Flin,X);

    % put matrices in the form f = A*y + B*u + G*p + d
    A = L(:,nu+1:nu+ny);
    B = L(:,1:nu);
    G = L(:,nu+ny+1:nx);
    Ai = L(:,nx+1:nx+ny);
    Af = L(:,nx+ny+1:nx+2*ny);
    d = -O;

    % vectorized matlab function
    E.A = sym2matrixfun(A,{t,PARAM},output);
    E.B = sym2matrixfun(B,{t,PARAM},output);
    E.G = sym2matrixfun(G,{t,PARAM},output);
    E.Ai = sym2matrixfun(Ai,{t,PARAM},output);
    E.Af = sym2matrixfun(Af,{t,PARAM},output);
    E.d = sym2matrixfun(d,{t,PARAM},output);

    % store indices (corresponding state equations)
    E.Ilin = find(Ilin);
    E.Inon = find(~Ilin);
else
    % treat all functions as nonlinear
    E.Ilin = [];
    E.Inon = 1:length(F);
end

% classify the functions as path or boundary types
if pathboundaryflag

    % sparsity pattern
    S = double(DF~=0);

    % continuous variable sparsity pattern
    Sc = S(:,1:(nu+ny));

    % time-varying dependency pattern
    St = hasSymType(DF,'symfunOf',t);

    % determine if a function is path or boundary
    E.pathboundary = any(Sc,2) | any(St,2);

else
    E.pathboundary = [];
end

% vectorized matlab function
E.f = sym2matrixfun(F,in1,output);

% vectorized matlab function
E.Df = sym2matrixfun(DF,in1,output);

% reshape into row vector
DFrow = reshape(DF.',1,[]);

% Jacobian again (second partial derivatives)
D2F = jacobian(DFrow,X');

% vectorized matlab function
D2f = sym2matrixfun(D2F,in1,output);

% reshape and store
if ~isempty(D2f)
    E.D2f = mat2cell(D2f,repmat(nx+2*ny,1,length(F)),nx+2*ny);
else
    E.D2f = [];
end

% store symbolic results
element.F = F;
element.DF = DF;
element.D2F = D2F;
E.element = element;

% (potentially) end the timer
if (opts.general.displevel > 0) % minimal
    opts.timer.sym = opts.timer.sym + toc(opts.timer.t2); % start timer
end

end

function [A,zeroflag] = sym2matrixfun(Ain,vars,output)

% return if input is empty
if isempty(Ain)
    A = [];
    zeroflag = true;
    return
end

% initialize with zeros
A = cell(size(Ain));
[A{:}] = deal(0);

% find nonzero entries
I = find(double(Ain~=0));

% return if all entries are zero
if isempty(I)
    A = [];
    zeroflag = true;
    return
else
    zeroflag = false;
end

% go through each nonzero entry
for k = I(:)'

    % extract
    g = Ain(k);

    % first see if it is a constant
    try
        A{k} = double(g);
    catch

        % construct appropriate matlab function
        if output == 1

            % convert to matlab function
            A{k} = matlabFunction(g,'Vars',vars);

        elseif output == 2

            % convert to matlab function
            Ai = matlabFunction(g,'Vars',vars);

            % interpolate optimization variables
            A{k} = @(t,T,X,param) Ai(t,param,interp1(T,X,t));

        end
    end
end

end