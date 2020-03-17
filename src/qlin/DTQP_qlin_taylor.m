%--------------------------------------------------------------------------
% DTQP_qlin_taylor.m
% Perform variable-order Taylor polynomial approximation on a matrix of
% functions with states, inputs, and parameters
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
function E = DTQP_qlin_taylor(f,form,in)

%--------------------------------------------------------------------------
% go through the options

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

% extra symbolic parameters
if isfield(in,'param')
	param = in.param;
else
    param = [];
end
%

if isfield(in,'output')
    output = in.output;
else
    output = 1; % vectorized matlab function
    % output = 2; % interp1-compatible matlab function
    % output = 3; % symbolic functions
end

%
switch form
    case 1
        norder = 2; % linear
        linflag = 0; % need linearization points?
    case 2
        norder = 3; % quadratic
        linflag = 0; % need linearization points?
    case 3
        norder = 2; % linear
        linflag = 1; % need linearization points?
    case 4
        norder = 3; % quadratic
        linflag = 1; % need linearization points?
    otherwise
        if isfield(in,'norder')
            norder = in.norder;
        else
            norder = in.norder;
        end
end

%--------------------------------------------------------------------------
% initialize the symbolic variable creation strings
strY = []; strU = []; strP = [];

% add states
for idx = 1:ny
    strY = [strY,'y',num2str(idx),' '];
end

% add controls
for idx = 1:nu
    strU = [strU,'u',num2str(idx),' '];
end

% add parameters
for idx = 1:np
    strP = [strP,'p',num2str(idx),' '];
end

% create symbolic variables
if linflag
   eval(['syms t dummy123 ',strU,upper(strU),strY,upper(strY),strP,upper(strP),param,' real'])
    eval(['AU = [',upper(strU),'];'])
    eval(['AY = [',upper(strY),'];'])
    eval(['AP = [',upper(strP),'];'])

else
    eval(['syms t dummy123 ',strU,strY,strP,param,' real'])
end

eval(['U = [',strU,'];'])
eval(['Y = [',strY,'];'])
eval(['P = [',strP,'];'])

if isempty(param)
    PARAM = dummy123;
else
	eval(['PARAM = [',param,'];'])
end

%--------------------------------------------------------------------------
% create the symbolic function
F = eval(f);

if linflag
    % Taylor polynomial approximation of F about X=A of order norder
    T = taylor(F,[U,Y,P],[AU,AY,AP],'Order',norder);
else
    % no Taylor approximation
    T = F;
end

% simplify the results
T = simplify(T);

switch form
    %----------------------------------------------------------------------
    case 1
        % convert linear equations to matrix notation
        [L,O] = equationsToMatrix(T,[Y,U,P]);

        % put matrices in the form f = A*y + B*u + G*p + d
        A = L(:,nu+1:nu+ny);
        B = L(:,1:nu);
        G = L(:,nu+ny+1:end);
        d = -O;

        % create the requested output structure
        if output == 1
            % vectorized matlab function
            E.A = sym2matrixfun(A,{t,PARAM},output);
            E.B = sym2matrixfun(B,{t,PARAM},output);
            E.G = sym2matrixfun(G,{t,PARAM},output);
            E.d = sym2matrixfun(d,{t,PARAM},output);
        elseif output == 3
            % symbolic functions
            E.A = A; E.B = B; E.G = G; E.d = d;
        end

    %----------------------------------------------------------------------
    case 2
        % put matrices in the form X'*h*X + g'*X + c
        h = hessian(T,[U,Y,P])/2; % hessian
        Tred = simplify(T-[U,Y,P]*h*[U,Y,P]'); % remove quadratic terms
        [g,c] = equationsToMatrix(Tred,[U,Y,P]); % gradient and constant

        % create the requested output structure
        if output == 1
            % vectorized matlab function
            E.H = sym2matrixfun(h,{t,PARAM},output);
            E.G = sym2matrixfun(g,{t,PARAM},output);
            E.C = sym2matrixfun(-c,{t,PARAM},output);
        elseif output == 3
            % symbolic functions
            E.H = h; E.G = g; E.C = -c;
        end

    %----------------------------------------------------------------------
    case 3
        % convert linear equations to matrix notation
        [L,O] = equationsToMatrix(T,[U,Y,P]);

        % put matrices in the form f = A*y + B*u + G*p + d
        A = L(:,nu+1:nu+ny);
        B = L(:,1:nu);
        G = L(:,nu+ny+1:end);
        d = -O;

        % create the requested output structure
        if output == 3
            % symbolic functions
            E.A = A; E.B = B; E.G = G; E.d = d;
        else % output == 1 or 2
            %
            E.A = sym2matrixfun(A,{t,PARAM,[AU,AY,AP]},output);
            E.B = sym2matrixfun(B,{t,PARAM,[AU,AY,AP]},output);
            E.G = sym2matrixfun(G,{t,PARAM,[AU,AY,AP]},output);
            E.d = sym2matrixfun(d,{t,PARAM,[AU,AY,AP]},output);
        end

    %----------------------------------------------------------------------
    case 4
        % put matrices in the form X'*h*X + g'*X + c
        h = hessian(T,[U,Y,P])/2; % hessian
        Tred = simplify(T-[U,Y,P]*h*[U,Y,P]','steps',1000); % remove quadratic terms
        [g,c] = equationsToMatrix(Tred,[U,Y,P]); % gradient and constant

        % create the requested output structure
        if output == 3
            % symbolic functions
            E.H = h; E.G = g; E.C = -c;
        else % output == 1 or 2
            % vectorized matlab function
            E.H = sym2matrixfun(h,{t,PARAM,[AY,AU,AP]},output);
            E.G = sym2matrixfun(g,{t,PARAM,[AY,AU,AP]},output);
            E.C = sym2matrixfun(-c,{t,PARAM,[AY,AU,AP]},output);
        end

    %----------------------------------------------------------------------
end

E.sym = T;

end

function A = sym2matrixfun(Ain,vars,output)

if isempty(Ain)
    A = [];
else
    % initialize
    A = cell(size(Ain));

    % go through each row
    for i = 1:size(Ain,1)
        % go through each column
        for j = 1:size(Ain,2)

            % get current function
            g = Ain(i,j);

            % first see if it is a constant
            try
                g = double(g);
                A{i,j} = g;
            catch
                if output == 1
                    Ai = matlabFunction(g,'Vars',vars);
                elseif output == 2
                    Ai = matlabFunction(g,'Vars',vars);
                    Ai = @(t,T,X,param) Ai(t,param,interp1(T,X,t));
                end
                A{i,j} = Ai;
            end
        end
    end
end

end